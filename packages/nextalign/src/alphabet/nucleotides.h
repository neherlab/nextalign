#pragma once

#include <nextalign/nextalign.h>

#include <string>

#include "../utils/to_underlying.h"

Nucleotide toNucleotide(char nuc);

char toChar(Nucleotide nuc);

inline std::ostream& operator<<(std::ostream& os, const Nucleotide& nucleotide) {
  os << std::string{to_underlying(nucleotide)};
  return os;
}
