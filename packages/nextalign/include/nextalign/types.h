#pragma once

#include <map>
#include <set>
#include <string>

struct NextalignOptions {
  std::set<std::string> genes;
};

struct Gene {
  std::string geneName;
  int start;
  int end;
  std::string strand;
  int frame;
};


using GeneMap = std::map<std::string, Gene>;


struct Alignment {
  std::string query;
  std::string ref;
  int alignmentScore;
};
