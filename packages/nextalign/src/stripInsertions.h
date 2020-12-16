#pragma once

#include <nextalign/types.h>

#include <string>
#include <vector>

struct StripInsertionsResult {
  std::string queryStripped;
  std::vector<Insertion> insertions;
};

StripInsertionsResult stripInsertions(const std::string& ref, const std::string& query);
