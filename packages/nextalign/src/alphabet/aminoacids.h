#pragma once

#include <string>

#include "../utils/to_underlying.h"
#include "letter.h"

constexpr const char CHAR_AMINOACID_UNKNOWN = 'X';
constexpr const char CHAR_AMINOACID_GAP = '-';
constexpr const char CHAR_AMINOACID_STOP = '*';

enum class Aminoacid : char {
  A = 0,
  B = 1,// D | N
  C = 2,
  D = 3,
  E = 4,
  F = 5,
  G = 6,
  H = 7,
  I = 8,
  J = 9,// L | I
  K = 10,
  L = 11,
  M = 12,
  N = 13,
  O = 14,// (rare) Pyrrolysine
  P = 15,
  Q = 16,
  R = 17,
  S = 18,
  T = 19,
  U = 20,// (rare) Selenocysteine
  V = 21,
  W = 22,
  Y = 23,
  Z = 24,// E | Q
  X = 25,
  STOP = 26,
  GAP = 27,
  SIZE,
};

using AminoacidSequence = Sequence<Aminoacid>;

Aminoacid toAminoacid(char aa);

char toChar(Aminoacid aa);

AminoacidSequence toAminoacidSequence(const std::string& seq);

std::string toString(const AminoacidSequence& seq);

inline std::ostream& operator<<(std::ostream& os, const Aminoacid& aminoacid) {
  os << std::string{to_underlying(aminoacid)};
  return os;
}
