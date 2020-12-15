#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

struct NextalignOptions {
  std::set<std::string> genes;
};

struct Gene {
  std::string geneName;
  int start;
  int end;
  std::string strand;
  int frame;
  int length;
};


using GeneMap = std::map<std::string, Gene>;


struct Alignment {
  std::string query;
  std::string ref;
  int alignmentScore;
};

struct Insertion {
  int begin;
  int end;
  std::string seq;
};

struct AlignmentImproved : public Alignment {
  std::vector<Insertion> insertions;
};
