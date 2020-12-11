#include "matchNuc.h"

#include <fmt/format.h>
#include <frozen/map.h>
#include <frozen/set.h>

#include <array>

namespace {
  class ErrorNucleotideInvalid : public std::runtime_error {
  public:
    explicit ErrorNucleotideInvalid(char nuc) : std::runtime_error(fmt::format("Invalid nucleotide: \"{:c}\"", nuc)) {}
  };


  static constexpr const frozen::set<IupacNucCodes, 4> canonicalNucleotides = /* NOLINT(cert-err58-cpp) */ {
    IupacNucCodes::A,
    IupacNucCodes::C,
    IupacNucCodes::G,
    IupacNucCodes::T,
  };


  // clang-format off
static constexpr const std::array<int, IupacNucCodes::SIZE*IupacNucCodes::SIZE> scoringMatrixNuc = {
  /*           00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15 */
  /*            U   T   A   W   C   Y   M   H   G   K   R   D   S   B   V   N */
  /* 00   U */  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 01   T */  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,  0,  1,
  /* 02   A */  0,  0,  1,  1,  0,  0,  1,  1,  0,  0,  1,  1,  0,  0,  1,  1,
  /* 03   W */  0,  1,  1,  1,  0,  1,  1,  1,  0,  1,  1,  1,  0,  1,  1,  1,
  /* 04   C */  0,  0,  0,  0,  1,  1,  1,  1,  0,  0,  0,  0,  1,  1,  1,  1,
  /* 05   Y */  0,  1,  0,  1,  1,  1,  1,  1,  0,  1,  0,  1,  1,  1,  1,  1,
  /* 06   M */  0,  0,  1,  1,  1,  1,  1,  1,  0,  0,  1,  1,  1,  1,  1,  1,
  /* 07   H */  0,  1,  1,  1,  1,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,
  /* 08   G */  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,
  /* 09   K */  0,  1,  0,  1,  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  /* 10   R */  0,  0,  1,  1,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  /* 11   D */  0,  1,  1,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  /* 12   S */  0,  0,  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  /* 13   B */  0,  1,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  /* 14   V */  0,  0,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  /* 15   N */  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  };
  // clang-format on


  static constexpr const frozen::map<char, IupacNucCodes, 16> nucToIupac = /* NOLINT(cert-err58-cpp) */ {
    /* 00 */ {'U', IupacNucCodes::U},
    /* 01 */ {'T', IupacNucCodes::T},
    /* 02 */ {'A', IupacNucCodes::A},
    /* 03 */ {'W', IupacNucCodes::W},
    /* 04 */ {'C', IupacNucCodes::C},
    /* 05 */ {'Y', IupacNucCodes::Y},
    /* 06 */ {'M', IupacNucCodes::M},
    /* 07 */ {'H', IupacNucCodes::H},
    /* 08 */ {'G', IupacNucCodes::G},
    /* 09 */ {'K', IupacNucCodes::K},
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

}// namespace

int lookupNucMatchScoreIupac(IupacNucCodes x, IupacNucCodes y) {
  return scoringMatrixNuc[x * IupacNucCodes::SIZE + y];//NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
}

int lookupNucMatchScore(char x, char y) {
  return lookupNucMatchScoreIupac(toIupac(x), toIupac(y));
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
  return lookupNucMatchScoreIupac(toIupac(x), toIupac(y)) > 0;
}
