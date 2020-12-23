#pragma once

#include <string>

#include "../utils/to_underlying.h"
#include "letter.h"

enum class Nucleotide : char {
  U = 0,
  T = 1,
  A = 2,
  W = 3,
  C = 4,
  Y = 5,
  M = 6,
  H = 7,
  G = 8,
  K = 9,
  R = 10,
  D = 11,
  S = 12,
  B = 13,
  V = 14,
  N = 15,
  GAP = 16,
  SIZE = 17,
};

using NucleotideSequence = Sequence<Nucleotide>;

Nucleotide toNucleotide(char nuc);

char toChar(Nucleotide nuc);

NucleotideSequence toNucleotideSequence(const std::string& seq);

std::string toString(const NucleotideSequence& seq);

inline std::ostream& operator<<(std::ostream& os, const Nucleotide& nucleotide) {
  os << std::string{to_underlying(nucleotide)};
  return os;
}
