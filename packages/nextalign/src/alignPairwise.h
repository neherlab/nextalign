#pragma once

#include <nextalign/types.h>

#include <string>
#include <vector>

#include "matchNuc.h"
#include "utils/vector2d.h"

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
  vector2d<int> scores;
  vector2d<int> paths;
};


struct AlignmentParameters {
  int gapExtend;
  int gapOpen;
  int gapClose;
  int misMatch;
  int match;
};


SeedAlignment seedAlignment(const std::string& query, const std::string& ref);

ForwardTrace scoreMatrix(const std::string& query, const std::string& ref, ScoreLookupFunction scoreLookupFunction, const std::vector<int>& gapOpenClose,
  int bandWidth, int meanShift);

Alignment backTrace(const std::string& query, const std::string& ref, const vector2d<int>& scores,
  const vector2d<int>& paths, int meanShift);

Alignment alignPairwise(
  const std::string& query, const std::string& ref, ScoreLookupFunction scoreLookupFunction, int minimalLength);

Alignment alignPairwise(
  const std::string& query, const std::string& ref, ScoreLookupFunction scoreLookupFunction, const std::vector<int>& gapOpenClose, int minimalLength);
