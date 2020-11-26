/**
 * Adopted with modifications from https://github.com/iosonofabio/seqanpy
 * Thanks Fabio Zanini @iosonofabio
 */

#include "seqan2.h"

#include <seqan2nowarn.h>

#include "helpers.h"


namespace {
  using SequenceIupac = seqan::String<seqan::Iupac>;
  using AlignIupac = seqan::Align<SequenceIupac, seqan::ArrayGaps>;

  using SequenceAminoacid = seqan::String<seqan::AminoAcid>;
  using AlignAminoacid = seqan::Align<SequenceAminoacid, seqan::ArrayGaps>;

  using ScoreSimple = seqan::Score<int, seqan::Simple>;
  using seqan::assignSource;
  using seqan::globalAlignment;
  using seqan::localAlignment;
  using seqan::resize;
  using seqan::row;
  using seqan::rows;

  template<typename Align>
  Align StringsToAlignment(const std::string& seq1, const std::string& seq2) {
    Align ali;
    resize(rows(ali), 2);
    assignSource(row(ali, 0), seq1);
    assignSource(row(ali, 1), seq2);
    return ali;
  }

  template<typename Align>
  std::tuple<std::string, std::string> AlignmentToStrings(const Align& ali) {
    std::stringstream ss;

    ss << row(ali, 0);
    const auto seq1 = ss.str();

    ss.clear();

    const auto seq2 = ss.str();
    ss << row(ali, 1);

    return std::make_tuple(seq1, seq2);
  }

  template<typename Align>
  AlignmentResult align_global(const std::string& seq1, const std::string& seq2, const AlignOptionsGlobal& options) {
    const auto& band = options.band;
    const auto& score_match = options.score_match;
    const auto& score_mismatch = options.score_mismatch;
    const auto& score_gapext = options.score_gapext;
    const auto& score_gapopen = options.score_gapopen;

    auto align = StringsToAlignment<Align>(seq1, seq2);

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

  template<typename Align>
  AlignmentResult align_overlap(const std::string& seq1, const std::string& seq2, const AlignOptionsOverlap& options) {
    const auto& band = options.band;
    const auto& score_match = options.score_match;
    const auto& score_mismatch = options.score_mismatch;
    const auto& score_gapext = options.score_gapext;
    const auto& score_gapopen = options.score_gapopen;
    const auto& cut_flanks = options.cut_flanks;

    // TODO(ivan-aksamentov): this parameter is not used, and is not present in other functions. Remove or find use?
    NA_UNUSED(cut_flanks);

    auto align = StringsToAlignment<Align>(seq1, seq2);

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

  template<typename Align>
  AlignmentResult align_ladder(const std::string& seq1, const std::string& seq2, const AlignOptionsLadder& options) {
    const auto& band = options.band;
    const auto& score_match = options.score_match;
    const auto& score_mismatch = options.score_mismatch;
    const auto& score_gapext = options.score_gapext;
    const auto& score_gapopen = options.score_gapopen;

    auto align = StringsToAlignment<Align>(seq1, seq2);

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

  template<typename Align>
  AlignmentResult align_local(const std::string& seq1, const std::string& seq2, const AlignOptionsLocal& options) {
    const auto& score_match = options.score_match;
    const auto& score_mismatch = options.score_mismatch;
    const auto& score_gapext = options.score_gapext;
    const auto& score_gapopen = options.score_gapopen;

    auto align = StringsToAlignment<Align>(seq1, seq2);

    // Align (best alignment only)
    const auto scoreSimple = ScoreSimple(score_match, score_mismatch, score_gapext, score_gapopen);
    int score = localAlignment(align, scoreSimple);

    const auto [ali1, ali2] = AlignmentToStrings(align);
    return {.score = score, .seq1 = ali1, .seq2 = ali2};
  }
}// namespace


// For nucleotides
[[maybe_unused]] AlignmentResult align_global(
  const std::string& seq1, const std::string& seq2, const AlignOptionsGlobal& options) {
  return align_global<AlignIupac>(seq1, seq2, options);
}

[[maybe_unused]] AlignmentResult align_overlap(
  const std::string& seq1, const std::string& seq2, const AlignOptionsOverlap& options) {
  return align_overlap<AlignIupac>(seq1, seq2, options);
}

[[maybe_unused]] AlignmentResult align_ladder(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLadder& options) {
  return align_ladder<AlignIupac>(seq1, seq2, options);
}

[[maybe_unused]] AlignmentResult align_local(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLocal& options) {
  return align_local<AlignIupac>(seq1, seq2, options);
}


// For aminoacids
[[maybe_unused]] AlignmentResult align_aa_global(
  const std::string& seq1, const std::string& seq2, const AlignOptionsGlobal& options) {
  return align_global<AlignAminoacid>(seq1, seq2, options);
}

[[maybe_unused]] AlignmentResult align_aa_overlap(
  const std::string& seq1, const std::string& seq2, const AlignOptionsOverlap& options) {
  return align_overlap<AlignAminoacid>(seq1, seq2, options);
}

[[maybe_unused]] AlignmentResult align_aa_ladder(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLadder& options) {
  return align_ladder<AlignAminoacid>(seq1, seq2, options);
}

[[maybe_unused]] AlignmentResult align_aa_local(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLocal& options) {
  return align_local<AlignAminoacid>(seq1, seq2, options);
}
