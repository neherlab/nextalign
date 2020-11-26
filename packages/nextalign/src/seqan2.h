/**
 * Adopted with modifications from https://github.com/iosonofabio/seqanpy
 * Thanks Fabio Zanini @iosonofabio
 */
#pragma once

#include <string>

struct AlignmentResult {
  int score;
  std::string seq1;
  std::string seq2;
};


struct AlignOptionsGlobal {
  int band;
  int score_match;
  int score_mismatch;
  int score_gapext;
  int score_gapopen;
};

struct AlignOptionsOverlap {
  int band;
  int score_match;
  int score_mismatch;
  int score_gapext;
  int score_gapopen;
  int cut_flanks;
};

struct AlignOptionsLadder {
  int band;
  int score_match;
  int score_mismatch;
  int score_gapext;
  int score_gapopen;
};

struct AlignOptionsLocal {
  int score_match;
  int score_mismatch;
  int score_gapext;
  int score_gapopen;
};


// For nucleotides
[[maybe_unused]] AlignmentResult align_global(
  const std::string& seq1, const std::string& seq2, const AlignOptionsGlobal& options);

[[maybe_unused]] AlignmentResult align_overlap(
  const std::string& seq1, const std::string& seq2, const AlignOptionsOverlap& options);

[[maybe_unused]] AlignmentResult align_ladder(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLadder& options);

[[maybe_unused]] AlignmentResult align_local(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLocal& options);


// For aminoacids
[[maybe_unused]] AlignmentResult align_aa_global(
  const std::string& seq1, const std::string& seq2, const AlignOptionsGlobal& options);

[[maybe_unused]] AlignmentResult align_aa_overlap(
  const std::string& seq1, const std::string& seq2, const AlignOptionsOverlap& options);

[[maybe_unused]] AlignmentResult align_aa_ladder(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLadder& options);

[[maybe_unused]] AlignmentResult align_aa_local(
  const std::string& seq1, const std::string& seq2, const AlignOptionsLocal& options);
