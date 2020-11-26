/**
 * Adopted with modifications from https://github.com/iosonofabio/seqanpy
 * Thanks Fabio Zanini @iosonofabio
 */

#include "seqan2.h"

#include <seqan2nowarn.h>

#include "helpers.h"


namespace {
  using Sequence = seqan::String<char>;
  using Align = seqan::Align<Sequence, seqan::ArrayGaps>;
  using ScoreSimple = seqan::Score<int, seqan::Simple>;
  using seqan::assignSource;
  using seqan::globalAlignment;
  using seqan::localAlignment;
  using seqan::resize;
  using seqan::row;
  using seqan::rows;

  Align StringsToAlignment(const std::string& seq1, const std::string& seq2) {
    Align ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), seq1);
    assignSource(row(ali, 1), seq2);
    return ali;
  }

  std::tuple<std::string, std::string> AlignmentToStrings(const Align& ali) {
    std::stringstream ss1, ss2;

    ss1 << row(ali, 0);
    const auto seq1 = ss1.str();

    ss2 << row(ali, 1);
    const auto seq2 = ss2.str();

    return std::make_tuple(seq1, seq2);
  }

}// namespace


[[maybe_unused]] AlignmentResult align_global(
  const std::string& seq1, const std::string& seq2, const AlignOptionsGlobal& options) {
  const auto& band = options.band;
  const auto& score_match = options.score_match;
  const auto& score_mismatch = options.score_mismatch;
  const auto& score_gapext = options.score_gapext;
  const auto& score_gapopen = options.score_gapopen;

  Align align = StringsToAlignment(seq1, seq2);

  const auto scoreSimple = ScoreSimple(score_match, score_mismatch, score_gapext, score_gapopen);
  int score;// NOLINT(cppcoreguidelines-init-variables)
  if (band >= 0) {
    score = globalAlignment(align, scoreSimple, -band, band);
  } else {
    score = globalAlignment(align, scoreSimple);
  }

  const auto [ali1, ali2] = AlignmentToStrings(align);
  return {.score = score, .seq1 = ali1, .seq2 = ali2};
}

[[maybe_unused]] AlignmentResult align_overlap(
  const std::string& seq1, const std::string& seq2, const AlignOptionsOverlap& options) {
  const auto& band = options.band;
  const auto& score_match = options.score_match;
  const auto& score_mismatch = options.score_mismatch;
  const auto& score_gapext = options.score_gapext;
  const auto& score_gapopen = options.score_gapopen;
  const auto& cut_flanks = options.cut_flanks;

  // TODO(ivan-aksamentov): this parameter is not used, and is not present in other functions. Remove or find use?
  NA_UNUSED(cut_flanks);

  Align align = StringsToAlignment(seq1, seq2);

  // Align
  const auto scoreSimple = ScoreSimple(score_match, score_mismatch, score_gapext, score_gapopen);
  using AlignConfigOverlap = seqan::AlignConfig<true, false, false, true>;
  int score;// NOLINT(cppcoreguidelines-init-variables)
  if (band >= 0) {
    score = globalAlignment(align, scoreSimple, AlignConfigOverlap(), -band, band);
  } else {
    score = globalAlignment(align, scoreSimple, AlignConfigOverlap());
  }

  const auto [ali1, ali2] = AlignmentToStrings(align);
  return {.score = score, .seq1 = ali1, .seq2 = ali2};
}

[[maybe_unused]] AlignmentResult align_ladder(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLadder& options) {
  const auto& band = options.band;
  const auto& score_match = options.score_match;
  const auto& score_mismatch = options.score_mismatch;
  const auto& score_gapext = options.score_gapext;
  const auto& score_gapopen = options.score_gapopen;

  Align align = StringsToAlignment(seq1, seq2);

  // Align
  const auto scoreSimple = ScoreSimple(score_match, score_mismatch, score_gapext, score_gapopen);
  using AlignConfigLadder = seqan::AlignConfig<true, false, true, false>;
  int score;// NOLINT(cppcoreguidelines-init-variables)
  if (band >= 0) {
    score = globalAlignment(align, scoreSimple, AlignConfigLadder(), -band, band);
  } else {
    score = globalAlignment(align, scoreSimple, AlignConfigLadder());
  }

  const auto [ali1, ali2] = AlignmentToStrings(align);
  return {.score = score, .seq1 = ali1, .seq2 = ali2};
}


[[maybe_unused]] AlignmentResult align_local(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLocal& options) {
  const auto& score_match = options.score_match;
  const auto& score_mismatch = options.score_mismatch;
  const auto& score_gapext = options.score_gapext;
  const auto& score_gapopen = options.score_gapopen;

  Align align = StringsToAlignment(seq1, seq2);

  // Align (best alignment only)
  const auto scoreSimple = ScoreSimple(score_match, score_mismatch, score_gapext, score_gapopen);
  int score = localAlignment(align, scoreSimple);

  const auto [ali1, ali2] = AlignmentToStrings(align);
  return {.score = score, .seq1 = ali1, .seq2 = ali2};
}
