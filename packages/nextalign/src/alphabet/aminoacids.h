#pragma once

#include <nextalign/nextalign.h>

#include <string>

#include "../utils/to_underlying.h"

Aminoacid toAminoacid(char aa);

char toChar(Aminoacid aa);

inline std::ostream& operator<<(std::ostream& os, const Aminoacid& aminoacid) {
  os << std::string{to_underlying(aminoacid)};
  return os;
}
