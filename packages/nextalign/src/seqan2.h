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

[[maybe_unused]] AlignmentResult align_global(const std::string& seq1, const std::string& seq2, const AlignOptionsGlobal& options);


struct AlignOptionsOverlap {
  int band;
  int score_match;
  int score_mismatch;
  int score_gapext;
  int score_gapopen;
  int cut_flanks;
};

AlignmentResult align_overlap(const std::string& seq1, const std::string& seq2, const AlignOptionsOverlap& options);

struct AlignOptionsLadder {
  int band;
  int score_match;
  int score_mismatch;
  int score_gapext;
  int score_gapopen;
};

AlignmentResult align_ladder(const std::string& seq1, const std::string& seq2, const AlignOptionsLadder& options);


struct AlignOptionsLocal {
  int score_match;
  int score_mismatch;
  int score_gapext;
  int score_gapopen;
};

[[maybe_unused]] AlignmentResult align_local(const std::string& seq1, const std::string& seq2, const AlignOptionsLocal& options);
