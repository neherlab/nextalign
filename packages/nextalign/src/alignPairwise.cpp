#include "alignPairwise.h"
#include <ctime>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <string>
#include <vector>
#include <iostream>
#include "matchNuc.h"

namespace details {
  int round(double x) {
    return static_cast<int>(std::round(x));
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


struct SeedMatch {
  int shift;
  int score;
};

struct SeedAlignment {
  int meanShift;
  int bandWidth;
};

struct ForwardTrace {
  std::vector<std::vector<int>> scores;
  std::vector<std::vector<int>> paths;
};


struct AlignmentParameters {
  int gapExtend;
  int gapOpen;
  int gapClose;
  int misMatch;
  int match;
};

AlignmentParameters alignmentParameters = {
  .gapExtend = 0,
  .gapOpen = -2,
  .gapClose = -2,
  .misMatch = -1,
  .match = 3,
};


// determine the position where a particular kmer (string of length k) matches the reference sequence
SeedMatch seedMatch(const std::string& kmer, const std::string& ref, const int start_pos) {
  int tmpScore = 0;
  int maxScore = 0;
  int maxShift = -1;
  for (int shift = start_pos; shift < ref.size() - kmer.size(); ++shift) {
    tmpScore = 0;
    for (int pos = 0; pos < kmer.size(); ++pos) {
      if (kmer[pos] == ref[shift + pos]) {
        tmpScore++;
      }
      // TODO: this speeds up seed-matching by disregarding bad seeds.
      if (tmpScore+3<pos){break;}
    }
    if (tmpScore > maxScore) {
      maxScore = tmpScore;
      maxShift = shift;
      // if maximal score is reached
      if (tmpScore==kmer.size()) {break;}
    }
  }
  return {.shift = maxShift, .score = maxScore};
}


SeedAlignment seedAlignment(const std::string& query, const std::string& ref) {
  constexpr const int nSeeds = 9;
  constexpr const int seedLength = 21;

  const int margin = ref.size() > 10000 ? 100 : details::round(ref.size() / 100.0);
  const int bandWidth = std::min(ref.size(), query.size());
  int start_pos = 0;
  if (bandWidth < 2 * seedLength) {
    return {.meanShift = 0, .bandWidth = bandWidth};
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
    const auto tmpMatch = seedMatch(query.substr(qPos, seedLength), ref, start_pos);

    // TODO: is this comment out of date?
    // only use seeds that match at least 90%
    if (tmpMatch.score >= 0.9 * seedLength) {
      seedMatches.push_back({qPos, tmpMatch.shift, tmpMatch.shift - qPos, tmpMatch.score});
      start_pos = tmpMatch.shift;
    }
    // std::cout <<qPos<<"  "<<tmpMatch.shift<< " "<<tmpMatch.score<<" "<<seedMatches[seedMatches.size()-1][2]<<"\n";
  }
  if (seedMatches.size() < 2) {
    throw ErrorAlignmentNoSeedMatches();
  }

  // given the seed matches, determine the maximal and minimal shifts
  // this shift is the typical amount the query needs shifting to match ref
  // ref:   ACTCTACTGC-TCAGAC
  // query: ----TCACTCATCT-ACACCGAT  => shift = 4, then 3, 4 again
  //TODO: This thing doesn't work:
  // const auto& [minShiftIt, maxShiftIt] = std::minmax(seedMatches.cbegin(), seedMatches.cend(),//
  //   [](const auto& left, const auto& right) { return left[2] < right[2]; });

  int minShift = ref.size();
  int maxShift = -ref.size();
  for (int si=0; si<seedMatches.size(); si++){
    if (seedMatches[si][2]<minShift){
      minShift = seedMatches[si][2];
    }
    if (seedMatches[si][2]>maxShift){
      maxShift = seedMatches[si][2];
    }
  }

  // std::cout <<minShift<<" "<<maxShift<<"\n";
  const int meanShift = details::round(0.5 * (minShift + maxShift));
  const int bandWidthFinal = maxShift - minShift + 9;
  return {.meanShift = meanShift, .bandWidth = bandWidthFinal};
}

// self made argmax function
template<typename Container>
std::pair<int, int> argmax(const Container& d) {
  int tmpmax = d[0];
  int tmpii = 0;

  for (int i = 1; i < d.size(); ++i) {
    const auto& x = d[i];
    if (x > tmpmax) {
      tmpmax = x;
      tmpii = i;
    }
  }

  return std::make_pair(tmpii, tmpmax);
}

ForwardTrace scoreMatrix(const std::string& query, const std::string& ref, int bandWidth, int meanShift) {
  // TODO: Avoid creating this lambda function
  const auto indexToShift = [&bandWidth, &meanShift]//
    (int si) {                                      //
      return si - bandWidth + meanShift;
    };

  // allocate a matrix to record the matches
  const int rowLength = ref.size() + 1;// TODO: Fix narrowing conversion
  std::vector<std::vector<int>> scores;// TODO: Avoid 2D vectors, use contiguous storage instead
  std::vector<std::vector<int>> paths;
  for (int shift = -bandWidth; shift < bandWidth + 1; ++shift) {
    scores.emplace_back(rowLength);
    paths.emplace_back(rowLength);
  }

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
  const int gapClose = alignmentParameters.gapClose;
  const int misMatch = alignmentParameters.misMatch;
  const int match = alignmentParameters.match;

  // TODO: Give these variables some meaningful names
  // TODO: Try to narrow the scope of these variables. Do all of these variables
  //  need to be forward-declared an uninitialized?
  const int END_OF_SEQUENCE = -1;
  int si;
  int ri;
  int shift;
  int tmpMatch;
  std::array<int, 4> cmp;
  int qPos;
  int origin;
  int score;
  int qGapOpen;
  int rGapOpen;
  for (ri = 0; ri < ref.size(); ri++) {
    for (si = 2 * bandWidth; si >= 0; si--) {
      shift = indexToShift(si);
      qPos = ri - shift;

      if (qPos < 0) {
        // preceeds query sequence -- no score, origin is query gap
        // we could fill all of this at once
        score = 0;
        origin = 3;
      } else if (qPos < query.size()) {
        // if the shifted position is within the query sequence
        tmpMatch = isMatch(query[qPos], ref[ri]) ? match : misMatch;
        // if the previous move included a gap, add gap-close cost
        if (paths[si][ri] == 2 || paths[si][ri] == 3) {
          tmpMatch += gapClose;
        }

        // determine whether the previous move was a reference or query gap
        rGapOpen = si < 2 * bandWidth ? (paths[si + 1][ri + 1] == 2 ? 0 : gapOpen) : 0;
        qGapOpen = si > 0 ? (paths[si - 1][ri] == 3 ? 0 : gapOpen) : 0;

        // calculate scores
        cmp = {
          0,                        // unaligned
          scores[si][ri] + tmpMatch,// match -- shift stays the same
          si < 2 * bandWidth ? scores[si + 1][ri + 1] + gapExtend + rGapOpen : gapExtend,// putting a gap into ref
          si > 0 ? scores[si - 1][ri] + gapExtend + qGapOpen : gapExtend,                // putting a gap into query
        };
        // determine best move and best score
        const auto [tmpOrigin, tmpScore] = argmax(cmp);
        origin = tmpOrigin;
        score = tmpScore;
      } else {
        // past query sequence -- mark as sequence end
        score = END_OF_SEQUENCE;
        origin = END_OF_SEQUENCE;
      }
      paths[si][ri + 1] = origin;
      scores[si][ri + 1] = score;
    }
  }

  return {.scores = scores, .paths = paths};
}

Alignment backTrace(const std::string& query, const std::string& ref, const std::vector<std::vector<int>>& scores,
  const std::vector<std::vector<int>>& paths, int meanShift) {
  const int bandWidth = (scores.size() - 1) / 2;
  const int rowLength = scores[0].size();

  // TODO: Avoid creating this lambda function
  const auto indexToShift = [&bandWidth, &meanShift]//
    (int si) {                                      //
      return si - bandWidth + meanShift;
    };


  std::vector<std::pair<char, char>> aln;
  std::string aln_ref, aln_query;
  aln_ref.reserve(rowLength + 3*bandWidth);
  aln_query.reserve(rowLength + 3*bandWidth);

  // Determine the best alignment by picking the optimal score at the end of the query
  // const lastIndexByShift = scores.map((d, i) = > Math.min(rowLength - 1, query.size() + indexToShift(i)));
  // const lastScoreByShift = scores.map((d, i) = > d[lastIndexByShift[i]]);


  std::vector<int> lastScoreByShift;
  std::vector<int> lastIndexByShift;
  lastScoreByShift.reserve(scores.size());
  lastIndexByShift.reserve(scores.size());

  int si = 0, bestScore = 0;
  for (int i=0; i<scores.size(); i++){
    lastIndexByShift[i] = rowLength - 1 < query.size() + indexToShift(i) ? rowLength - 1 : query.size() + indexToShift(i);
    lastScoreByShift[i] = scores[i][lastIndexByShift[i]];
    if (lastScoreByShift[i]>bestScore){
      bestScore = lastScoreByShift[i];
      si = i;
    }
  }

  const int shift = indexToShift(si);
  int origin;

  // determine position tuple qPos, rPos corresponding to the place it the matrix
  int rPos = lastIndexByShift[si] - 1;
  int qPos = rPos - shift;
  // add right overhang, i.e. unaligned parts of the query or reference the right end
  if (rPos < ref.size() - 1) {
    for (int ii = ref.size() - 1; ii > rPos; ii--) {
      aln_query += '-';
      aln_ref += ref[ii];
    }
  } else if (qPos < query.size() - 1) {
    for (int ii = query.size() - 1; ii > qPos; ii--) {
      aln_query += query[ii];
      aln_ref += '-';
    }
  }

  // do backtrace for aligned region
  while (rPos > 0 && qPos > 0) {
    // tmpMatch = ref[rPos] === query[qPos] || query[qPos] === "N" ? match : misMatch;
    origin = paths[si][rPos + 1];
    if (origin == 1) {
      // match -- decrement both strands and add match to alignment
      aln_query += query[qPos];
      aln_ref += ref[rPos];
      qPos--;
      rPos--;
    } else if (origin == 2) {
      // insertion in query -- decrement query, increase shift
      aln_query += query[qPos];
      aln_ref += '-';
      qPos--;
      si++;
    } else if (origin == 3) {
      // deletion in query -- decrement reference, reduce shift
      aln_query += '-';
      aln_ref += ref[rPos];
      rPos--;
      si--;
    } else {
      break;
    }
  }
  // add the last match
  aln_query += query[qPos];
  aln_ref += ref[rPos];

  // add left overhang
  if (rPos > 0) {
    for (int ii = rPos - 1; ii >= 0; ii--) {
      aln_query += '-';
      aln_ref += ref[ii];
    }
  } else if (qPos > 0) {
    for (int ii = qPos - 1; ii >= 0; ii--) {
      aln_query += query[ii];
      aln_ref += '-';
    }
  }

  // reverse and make sequence
  // std::reverse(aln.begin(), aln.end());
  std::reverse(aln_query.begin(), aln_query.end());
  std::reverse(aln_ref.begin(), aln_ref.end());

  // TODO: these caused errors for me -- eliminated the pair vector
  // const auto queryFinal = std::string(aln.size(), '-');
  // std::accumulate(aln.cbegin(), aln.cend(), query.begin(),//
  //   [](const auto& x, std::string& res) { return res + x[0]; });

  // const auto refFinal = std::string(aln.size(), '-');
  // std::accumulate(aln.cbegin(), aln.cend(), query.begin(),//
  //   [](const auto& x, std::string& res) { return res + x[1]; });

  return {
    .query = aln_query,
    .ref = aln_ref,
    .alignmentScore = bestScore,
  };
}

Alignment alignPairwise(const std::string& query, const std::string& ref, int minimalLength) {
  if (query.size() < minimalLength) {
    throw ErrorAlignmentSequenceTooShort();
  }
  clock_t t1, t2;

  // perform a number of seed matches to determine te rough alignment of query rel to ref
  t1 = std::clock();
  const SeedAlignment& seedAlignmentResult = seedAlignment(query, ref);
  const auto& bandWidth = seedAlignmentResult.bandWidth;
  const auto& meanShift = seedAlignmentResult.meanShift;
  // std::cout <<"shift "<<meanShift<<" band "<<bandWidth<<"\n";
  if (bandWidth > 400) {
    throw ErrorAlignmentBadSeedMatches();
  }

  t2 = std::clock();
  std::cout<<"\nseed matching: "<<t2-t1<<"\n";
  const ForwardTrace& forwardTrace = scoreMatrix(query, ref, bandWidth, meanShift);
  t1 = std::clock();
  std::cout<<"forward trace: "<<t1-t2<<"\n";
  const auto& scores = forwardTrace.scores;
  const auto& paths = forwardTrace.paths;

  return backTrace(query, ref, scores, paths, meanShift);
}
