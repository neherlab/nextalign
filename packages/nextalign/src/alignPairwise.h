#pragma once

#include <nextalign/types.h>

#include <string>
#include <vector>

struct NextalignOptions;

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


SeedAlignment seedAlignment(const std::string& query, const std::string& ref);

ForwardTrace scoreMatrix(const std::string& query, const std::string& ref, int bandWidth, int meanShift);

Alignment backTrace(const std::string& query, const std::string& ref, const std::vector<std::vector<int>>& scores,
  const std::vector<std::vector<int>>& paths, int meanShift);

Alignment alignPairwise(const std::string& query, const std::string& ref, const int min_length);
