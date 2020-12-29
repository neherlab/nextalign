#pragma once

#include <nextalign/nextalign.h>

#include <gsl/string_span>
#include <string>

#include "../nextalign_private.h"
#include "../utils/to_underlying.h"

using AminoacidSequenceSpan = SequenceSpan<Aminoacid>;

Aminoacid toAminoacid(char aa);

char toChar(Aminoacid aa);

inline std::ostream& operator<<(std::ostream& os, const Aminoacid& aminoacid) {
  os << std::string{to_underlying(aminoacid)};
  return os;
}
