#include "alignPairwise.h"

#include <algorithm>
#include <cmath>
#include <ctime>
#include <gsl/gsl>
#include <iostream>
#include <string>
#include <vector>

#include "../alphabet/aminoacids.h"
#include "../alphabet/nucleotides.h"
#include "../match/matchAa.h"
#include "../match/matchNuc.h"
#include "../utils/safe_cast.h"
#include "../utils/vector2d.h"


namespace details {
  int round(double x) {
    return safe_cast<int>(std::round(x));
  }
}// namespace details


class ErrorAlignmentNoSeedMatches : public std::runtime_error {
public:
  explicit ErrorAlignmentNoSeedMatches() : std::runtime_error("Unable to align: no seed matches") {}
};

class ErrorAlignmentSequenceTooShort : public std::runtime_error {
public:
  explicit ErrorAlignmentSequenceTooShort() : std::runtime_error("Unable to align: sequence is too short") {}
};

class ErrorAlignmentBadSeedMatches : public std::runtime_error {
public:
  explicit ErrorAlignmentBadSeedMatches()
      : std::runtime_error("Unable to align: too many insertions, deletions, duplications, or ambiguous seed matches") {
  }
};


AlignmentParameters alignmentParameters = {
  .gapExtend = 0,
  .gapOpen = -6,
  .misMatch = -1,
  .match = 3,
};

// store direction info for backtrace as bits in paths matrix
// these indicate the currently optimal move
int MATCH=1<<0;
int refGAPmatrix=1<<1;
int qryGAPmatrix=1<<2;
// these are the override flags for gap extension
int refGAPextend=1<<3;
int qryGAPextend=1<<4;


// determine the position where a particular kmer (string of length k) matches the reference sequence
template<typename Letter>
SeedMatch seedMatch(
  const Sequence<Letter>& kmer, const Sequence<Letter>& ref, const int start_pos, const int allowed_mismatches) {
  const int refSize = safe_cast<int>(ref.size());
  const int kmerSize = safe_cast<int>(kmer.size());
  int tmpScore = 0;
  int maxScore = 0;
  int maxShift = -1;
  for (int shift = start_pos; shift < refSize - kmerSize; ++shift) {
    tmpScore = 0;
    for (int pos = 0; pos < kmerSize; ++pos) {
      if (kmer[pos] == ref[shift + pos]) {
        tmpScore++;
      }
      // TODO: this speeds up seed-matching by disregarding bad seeds.
      if (tmpScore + allowed_mismatches < pos) {
        break;
      }
    }
    if (tmpScore > maxScore) {
      maxScore = tmpScore;
      maxShift = shift;
      // if maximal score is reached
      if (tmpScore == kmerSize) {
        break;
      }
    }
  }
  return {.shift = maxShift, .score = maxScore};
}

template<typename Letter>
SeedAlignment seedAlignment(const Sequence<Letter>& query, const Sequence<Letter>& ref) {
  constexpr const int nSeeds = 9;
  constexpr const int seedLength = 21;
  constexpr const int allowed_mismatches = 3;

  const int margin = ref.size() > 10000 ? 100 : details::round(ref.size() / 100.0);
  const int bandWidth = details::round((ref.size()+query.size())*0.5) - 3;
  int start_pos = 0;
  if (bandWidth < 2 * seedLength) {
    return {.meanShift = details::round(((int)ref.size()-(int)query.size())*0.5), .bandWidth = bandWidth};
  }

  // TODO; give a name to this type.
  //  Maybe use something other than array? A struct with named fields to make
  //  the code in the end of the function less confusing?
  using TodoNameTheType = std::array<int, 4>;
  std::vector<TodoNameTheType> seedMatches;
  for (int ni = 0; ni < nSeeds; ++ni) {

    // TODO: give this variable a name
    // generate kmers equally spaced on the query
    const auto todoGiveAName = static_cast<double>(query.size() - seedLength - 2 * margin);
    const int qPos = details::round(margin + ((todoGiveAName) / (nSeeds - 1)) * ni);

    // TODO: verify that the `query.substr()` behavior is the same as JS `string.substr()`
    const auto tmpMatch = seedMatch(query.substr(qPos, seedLength), ref, start_pos, allowed_mismatches);

    // only use seeds with at most allowed_mismatches
    if (tmpMatch.score >= seedLength - allowed_mismatches) {
      seedMatches.push_back({qPos, tmpMatch.shift, tmpMatch.shift - qPos, tmpMatch.score});
      start_pos = tmpMatch.shift;
    }
  }
  if (seedMatches.size() < 2) {
    throw ErrorAlignmentNoSeedMatches();
  }

  // given the seed matches, determine the maximal and minimal shifts
  // this shift is the typical amount the query needs shifting to match ref
  // ref:   ACTCTACTGC-TCAGAC
  // query: ----TCACTCATCT-ACACCGAT  => shift = 4, then 3, 4 again

  int minShift = ref.size();
  int maxShift = -ref.size();
  for (auto& seedMatch : seedMatches) {
    if (seedMatch[2] < minShift) {
      minShift = seedMatch[2];
    }
    if (seedMatch[2] > maxShift) {
      maxShift = seedMatch[2];
    }
  }

  const int meanShift = details::round(0.5 * (minShift + maxShift));
  const int bandWidthFinal = maxShift - minShift + 9;
  return {.meanShift = meanShift, .bandWidth = bandWidthFinal};
}

template<typename Letter>
ForwardTrace scoreMatrix(const Sequence<Letter>& query, const Sequence<Letter>& ref, int bandWidth,
  int meanShift) {                                  // TODO: Avoid creating this lambda function
  const auto indexToShift = [&bandWidth, &meanShift]//
    (int si) {                                      //
      return si - bandWidth + meanShift;
    };
  // allocate a matrix to record the matches
  const int querySize = safe_cast<int>(query.size());
  const int refSize = safe_cast<int>(ref.size());
  const int n_rows = bandWidth * 2 + 1;
  const int n_cols = refSize + 1;

  vector2d<int> paths(n_rows, n_cols);
  // these could be reduced to vectors
  vector2d<int> scores(n_rows, n_cols);
  vector2d<int> refGaps(n_rows, n_cols);
  vector2d<int> qryGaps(n_rows, n_cols);

  // fill scores with alignment scores
  // The inner index scores[][ri] is the index of the reference sequence
  // the outer index si index the shift, together they define rPos=ri and qPos = ri-shift
  // if the colon marks the position in the sequence before rPos,qPos
  // R: ...ACT:X
  // Q: ...ACT:Y
  // 1) if X and Y are bases they either match or mismatch. shift doesn't change, rPos and qPos advance
  //    -> right horizontal step in the matrix
  // 2) if X is '-' and Y is a base, rPos stays the same and the shift decreases
  //    -> vertical step in the matrix from si+1 to si
  // 2) if X is a base and Y is '-', rPos advances the same and the shift increases
  //    -> diagonal step in the matrix from (ri,si-1) to (ri+1,si)
  const int gapExtend = alignmentParameters.gapExtend;
  const int gapOpen = alignmentParameters.gapOpen;
  const int misMatch = alignmentParameters.misMatch;
  const int match = alignmentParameters.match;

  // TODO: Give these variables some meaningful names
  // TODO: Try to narrow the scope of these variables. Do all of these variables
  //  need to be forward-declared an uninitialized?
  const int END_OF_SEQUENCE = -1;

#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-init-variables"
#pragma ide diagnostic ignored "cppcoreguidelines-pro-type-member-init"
  int si;
  int ri;
  int shift;
  int tmpMatch;
  int qPos;
  int origin;
  bool originRefGap, originQryGap;
  int score, tmpScore;
  int qGapOpen, qGapExtend;
  int rGapOpen, rGapExtend;
#pragma clang diagnostic pop
  for (si = 2 * bandWidth; si > bandWidth; si--) {
    paths(si,0) = qryGAPmatrix;
  }
  paths(bandWidth,0) = MATCH;
  for (si = bandWidth-1; si >=0; si--) {
    paths(si,0) = refGAPmatrix;
    qryGaps(si,0) = gapOpen;
  }
  for (ri = 0; ri < refSize; ri++) {
    for (si = 2 * bandWidth; si >= 0; si--) {
      shift = indexToShift(si);
      qPos = ri - shift;
      if (qPos < 0) {
        // precedes query sequence -- no score, origin is query gap
        // we could fill all of this at once
        score = 0;
        paths(si, ri+1) += qryGAPextend;
        refGaps(si,ri+1) = gapOpen;
        origin = qryGAPmatrix;
      } else if (qPos < querySize) {
        // if the shifted position is within the query sequence
        tmpMatch = lookupMatchScore(query[qPos], ref[ri]) > 0 ? match : misMatch;

        score = scores(si, ri) + tmpMatch;
        origin = MATCH;
        // determine whether the previous move was a reference or query gap
        if (si < 2 * bandWidth){
          rGapExtend = refGaps(si + 1, ri + 1) + gapExtend;
          rGapOpen = scores(si + 1, ri + 1) + gapOpen;
          if (rGapExtend>rGapOpen){
            tmpScore = rGapExtend;
            paths(si, ri+1) += refGAPextend;
          }else{
            tmpScore = rGapOpen;
          }
          refGaps(si, ri+1) = tmpScore;
          if (score < tmpScore){
            score = tmpScore;
            origin = refGAPmatrix;
          }
        }

        if (si > 0){
          qGapExtend = qryGaps(si - 1, ri) + gapExtend;
          qGapOpen = scores(si - 1, ri) + gapOpen;
          tmpScore = qGapExtend>qGapOpen ? qGapExtend : qGapOpen;
          if (qGapExtend>qGapOpen){
            tmpScore = qGapExtend;
            paths(si, ri+1) += qryGAPextend;
          }else{
            tmpScore = qGapOpen;
          }
          qryGaps(si, ri+1) = tmpScore;
          if (score < tmpScore){
            score = tmpScore;
            origin = qryGAPmatrix;
          }
        }
      } else {
        // past query sequence -- mark as sequence end
        score = END_OF_SEQUENCE;
        origin = END_OF_SEQUENCE;
      }
      // std::cout <<origin<<":"<<(paths(si, ri + 1)&refGAPextend)<<":"<<(paths(si, ri + 1)&qryGAPextend)<<":"<<score<<"\t";
      paths(si, ri + 1) += origin;
      scores(si, ri + 1) = score;
    }
    // std::cout <<"\n";
  }

  return {.scores = scores, .paths = paths};
}

template<typename Letter>
AlignmentResult<Letter> backTrace(const Sequence<Letter>& query, const Sequence<Letter>& ref,
  const vector2d<int>& scores, const vector2d<int>& paths, int meanShift) {
  const int rowLength = scores.num_cols();
  const int scoresSize = safe_cast<int>(scores.num_rows());
  const int querySize = safe_cast<int>(query.size());
  const int refSize = safe_cast<int>(ref.size());
  const int bandWidth = (scoresSize - 1) / 2;
  int currentMatrix = 0;

  // TODO: Avoid creating this lambda function
  const auto indexToShift = [&bandWidth, &meanShift]//
    (int si) {                                      //
      return si - bandWidth + meanShift;
    };


  std::vector<std::pair<char, char>> aln;
  Sequence<Letter> aln_ref;
  Sequence<Letter> aln_query;
  aln_ref.reserve(rowLength + 3 * bandWidth);
  aln_query.reserve(rowLength + 3 * bandWidth);

  // const lastIndexByShift = scores.map((d, i) = > Math.min(rowLength - 1, query.size() + indexToShift(i)));
  // const lastScoreByShift = scores.map((d, i) = > d[lastIndexByShift[i]]);


  std::vector<int> lastScoreByShift;
  std::vector<int> lastIndexByShift;
  lastScoreByShift.resize(scores.num_rows());
  lastIndexByShift.resize(scores.num_rows());

  // Determine the best alignment by picking the optimal score at the end of the query
  int si = 0;
  int bestScore = 0;
  for (int i = 0; i < scoresSize; i++) {
    lastIndexByShift[i] = rowLength - 1 < querySize + indexToShift(i) ? rowLength - 1 : query.size() + indexToShift(i);
    lastScoreByShift[i] = scores(i, lastIndexByShift[i]);
    if (lastScoreByShift[i] > bestScore) {
      bestScore = lastScoreByShift[i];
      si = i;
    }
  }

  const int shift = indexToShift(si);
  int origin;//NOLINT(cppcoreguidelines-init-variables)

  // determine position tuple qPos, rPos corresponding to the place it the matrix
  int rPos = lastIndexByShift[si] - 1;
  int qPos = rPos - shift;
  // add right overhang, i.e. unaligned parts of the query or reference the right end
  if (rPos < refSize - 1) {
    for (int ii = refSize - 1; ii > rPos; ii--) {
      aln_query += Letter::GAP;
      aln_ref += ref[ii];
    }
  } else if (qPos < querySize - 1) {
    for (int ii = querySize - 1; ii > qPos; ii--) {
      aln_query += query[ii];
      aln_ref += Letter::GAP;
    }
  }

  // do backtrace for aligned region
  while (rPos >= 0 && qPos >= 0) {
    origin = paths(si, rPos + 1);
    // std::cout<<si<<" "<<rPos<<" "<<origin<<" "<<currentMatrix<<"\n";
    if (origin&MATCH && currentMatrix==0) {
      // match -- decrement both strands and add match to alignment
      aln_query += query[qPos];
      aln_ref += ref[rPos];
      qPos--;
      rPos--;
    } else if ((origin&refGAPmatrix && currentMatrix==0) || currentMatrix==refGAPmatrix) {
      // insertion in ref -- decrement query, increase shift
      aln_query += query[qPos];
      aln_ref += Letter::GAP;
      qPos--;
      si++;
      if (origin&refGAPextend){
        // remain in gap-extension mode and ignore best-overall score
        currentMatrix=refGAPmatrix;
      }else{
        // close gap, return to best-overall score
        currentMatrix = 0;
      }
    } else if ((origin&qryGAPmatrix && currentMatrix==0) || currentMatrix==qryGAPmatrix) {
      // deletion in query -- decrement reference, reduce shift
      aln_query += Letter::GAP;
      aln_ref += ref[rPos];
      rPos--;
      si--;
      if (origin&qryGAPextend){
        // remain in gap-extension mode and ignore best-overall score
        currentMatrix=qryGAPmatrix;
      }else{
        // close gap, return to best-overall score
        currentMatrix = 0;
      }
    } else {
      break;
    }
  }

  // add left overhang
  if (rPos >= 0) {
    for (int ii = rPos; ii >= 0; ii--) {
      aln_query += Letter::GAP;
      aln_ref += ref[ii];
    }
  } else if (qPos >= 0) {
    for (int ii = qPos; ii >= 0; ii--) {
      aln_query += query[ii];
      aln_ref += Letter::GAP;
    }
  }

  // reverse and make sequence
  // std::reverse(aln.begin(), aln.end());
  std::reverse(aln_query.begin(), aln_query.end());
  std::reverse(aln_ref.begin(), aln_ref.end());

  return {
    .query = aln_query,
    .ref = aln_ref,
    .alignmentScore = bestScore,
  };
}

struct AlignPairwiseTag {};

template<typename Letter>
AlignmentResult<Letter> alignPairwise(
  const Sequence<Letter>& query, const Sequence<Letter>& ref, int minimalLength, AlignPairwiseTag) {
  const int querySize = query.size();
  if (querySize < minimalLength) {
    throw ErrorAlignmentSequenceTooShort();
  }

  const SeedAlignment& seedAlignmentResult = seedAlignment(query, ref);
  const auto& bandWidth = seedAlignmentResult.bandWidth;
  const auto& meanShift = seedAlignmentResult.meanShift;

  if (bandWidth > 400) {
    throw ErrorAlignmentBadSeedMatches();
  }
  const ForwardTrace& forwardTrace = scoreMatrix(query, ref, bandWidth, meanShift);
  const auto& scores = forwardTrace.scores;
  const auto& paths = forwardTrace.paths;

  return backTrace(query, ref, scores, paths, meanShift);
}


NucleotideAlignmentResult alignPairwise(
  const NucleotideSequence& query, const NucleotideSequence& ref, int minimalLength) {
  return alignPairwise(query, ref, minimalLength, AlignPairwiseTag{});
}

AminoacidAlignmentResult alignPairwise(
  const AminoacidSequence& query, const AminoacidSequence& ref, int minimalLength) {
  return alignPairwise(query, ref, minimalLength, AlignPairwiseTag{});
}
