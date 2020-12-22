#pragma once

#include <nextalign/types.h>

#include <string>
#include <vector>

#include "alphabet/nucleotides.h"

struct InsertionInternal {
  int begin;
  int end;
  NucleotideSequence seq;
};

struct StripInsertionsResult {
  NucleotideSequence queryStripped;
  std::vector<InsertionInternal> insertions;
};

StripInsertionsResult stripInsertions(const NucleotideSequence& ref, const NucleotideSequence& query);
