#include "matchNuc.h"

#include <fmt/format.h>

#include <array>
#include <map>
#include <set>

class ErrorNucleotideInvalid : public std::runtime_error {
public:
  explicit ErrorNucleotideInvalid(char nuc) : std::runtime_error(fmt::format("Invalid nucleotide: \"{:c}\"", nuc)) {}
};


const std::set<IupacNucCodes> canonicalNucleotides = /* NOLINT(cert-err58-cpp) */ {
  IupacNucCodes::A,
  IupacNucCodes::C,
  IupacNucCodes::G,
  IupacNucCodes::T,
};


// clang-format off
constexpr const std::array<int, IupacNucCodes::SIZE*IupacNucCodes::SIZE> scoringMatrixNuc = {
  /*       U   T   A   W   C   Y   M   H   G   K   R   D   S   B   V   N */
  /* U */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* T */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* A */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* W */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* C */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* Y */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* M */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* H */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* G */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* K */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* R */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* D */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* S */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* B */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* V */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
  /* N */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
};
// clang-format on


const std::map<char, IupacNucCodes> nucToIupac = /* NOLINT(cert-err58-cpp) */ {
  /* 0 */ {'U', IupacNucCodes::U},
  /* 1 */ {'T', IupacNucCodes::T},
  /* 2 */ {'A', IupacNucCodes::A},
  /* 3 */ {'W', IupacNucCodes::W},
  /* 4 */ {'C', IupacNucCodes::C},
  /* 5 */ {'Y', IupacNucCodes::Y},
  /* 6 */ {'M', IupacNucCodes::M},
  /* 7 */ {'H', IupacNucCodes::H},
  /* 8 */ {'G', IupacNucCodes::G},
  /* 9 */ {'K', IupacNucCodes::K},
  /* 10 */ {'R', IupacNucCodes::R},
  /* 11 */ {'D', IupacNucCodes::D},
  /* 12 */ {'S', IupacNucCodes::S},
  /* 13 */ {'B', IupacNucCodes::B},
  /* 14 */ {'V', IupacNucCodes::V},
  /* 15 */ {'N', IupacNucCodes::N},
};


IupacNucCodes toIupac(char nuc) {
  const auto it = nucToIupac.find(nuc);
  if (it == nucToIupac.end()) {
    throw ErrorNucleotideInvalid(nuc);
  }
  return it->second;
}


int lookupNucMatchScore(IupacNucCodes x, IupacNucCodes y) {
  return scoringMatrixNuc[x * IupacNucCodes::SIZE + y];//NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
}

bool isMatch(char x, char y) {
  if (x == y) {
    return true;
  }

  if (y == 'N') {
    return true;
  }

  if (x == 'N') {
    return true;
  }
  return lookupNucMatchScore(toIupac(x), toIupac(y)) > 0;
}
