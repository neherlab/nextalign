#pragma once

#include <map>
#include <set>
#include <string>


struct NextalignOptions {
  std::set<std::string> genes;
};

struct Gene {};

using GeneMap = std::map<std::string, Gene>;
