#pragma once

#include <functional>

enum IupacNucCodes {
  U [[maybe_unused]] = 0,
  T [[maybe_unused]] = 1,
  A [[maybe_unused]] = 2,
  W [[maybe_unused]] = 3,
  C [[maybe_unused]] = 4,
  Y [[maybe_unused]] = 5,
  M [[maybe_unused]] = 6,
  H [[maybe_unused]] = 7,
  G [[maybe_unused]] = 8,
  K [[maybe_unused]] = 9,
  R [[maybe_unused]] = 10,
  D [[maybe_unused]] = 11,
  S [[maybe_unused]] = 12,
  B [[maybe_unused]] = 13,
  V [[maybe_unused]] = 14,
  N [[maybe_unused]] = 15,
  SIZE = 16,
};

using ScoreLookupFunction = std::function<int(char, char)>;

int lookupNucMatchScore(char x, char y);
