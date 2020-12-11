#include "matchAa.h"

#include <fmt/format.h>
#include <frozen/map.h>
#include <frozen/set.h>

#include <array>

#include "decode.h"

namespace {
  class ErrorAminoacidInvalid : public std::runtime_error {
  public:
    explicit ErrorAminoacidInvalid(char aa) : std::runtime_error(fmt::format("Invalid aminoacid: \"{:c}\"", aa)) {}
  };


  enum AminoacidCode {
    A [[maybe_unused]] = 0,
    B [[maybe_unused]] = 1,// D | N
    C [[maybe_unused]] = 2,
    D [[maybe_unused]] = 3,
    E [[maybe_unused]] = 4,
    F [[maybe_unused]] = 5,
    G [[maybe_unused]] = 6,
    H [[maybe_unused]] = 7,
    I [[maybe_unused]] = 8,
    J [[maybe_unused]] = 9,// L | I
    K [[maybe_unused]] = 10,
    L [[maybe_unused]] = 11,
    M [[maybe_unused]] = 12,
    N [[maybe_unused]] = 13,
    O [[maybe_unused]] = 14,// (rare) Pyrrolysine
    P [[maybe_unused]] = 15,
    Q [[maybe_unused]] = 16,
    R [[maybe_unused]] = 17,
    S [[maybe_unused]] = 18,
    T [[maybe_unused]] = 19,
    U [[maybe_unused]] = 20,// (rare) Selenocysteine
    V [[maybe_unused]] = 21,
    W [[maybe_unused]] = 22,
    Y [[maybe_unused]] = 23,
    Z [[maybe_unused]] = 24,// E | Q
    X [[maybe_unused]] = 25,
    STOP [[maybe_unused]] = 26,
    SIZE,
  };

  static constexpr const frozen::map<Aminoacid, AminoacidCode, 27> charToAminoacid = /* NOLINT(cert-err58-cpp) */ {
    /* 00 */ {'A', AminoacidCode::A},
    /* 01 */ {'B', AminoacidCode::B},
    /* 02 */ {'C', AminoacidCode::C},
    /* 03 */ {'D', AminoacidCode::D},
    /* 04 */ {'E', AminoacidCode::E},
    /* 05 */ {'F', AminoacidCode::F},
    /* 06 */ {'G', AminoacidCode::G},
    /* 07 */ {'H', AminoacidCode::H},
    /* 08 */ {'I', AminoacidCode::I},
    /* 09 */ {'J', AminoacidCode::J},
    /* 10 */ {'K', AminoacidCode::K},
    /* 11 */ {'L', AminoacidCode::L},
    /* 12 */ {'M', AminoacidCode::M},
    /* 13 */ {'N', AminoacidCode::N},
    /* 14 */ {'O', AminoacidCode::O},
    /* 15 */ {'P', AminoacidCode::P},
    /* 16 */ {'Q', AminoacidCode::Q},
    /* 17 */ {'R', AminoacidCode::R},
    /* 18 */ {'S', AminoacidCode::S},
    /* 19 */ {'T', AminoacidCode::T},
    /* 20 */ {'U', AminoacidCode::U},
    /* 21 */ {'V', AminoacidCode::V},
    /* 22 */ {'W', AminoacidCode::W},
    /* 23 */ {'Y', AminoacidCode::Y},
    /* 24 */ {'Z', AminoacidCode::Z},
    /* 25 */ {AMINOACID_UNKNOWN, AminoacidCode::X},
    /* 26 */ {AMINOACID_STOP, AminoacidCode::STOP},
  };


  // clang-format off
static constexpr const std::array<int, AminoacidCode::SIZE*AminoacidCode::SIZE> scoringMatrixAa = {
  /*           00  01  02  03  04  05  06  07  08  09  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26 */
  /*            A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   Y   Z   *   X */
  /* 00   A */  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 01   B */  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 02   C */  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 03   D */  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 04   E */  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,
  /* 05   F */  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 06   G */  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 07   H */  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 08   I */  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 09   J */  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 10   K */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 11   L */  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 12   M */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 13   N */  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 14   O */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 15   P */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 16   Q */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,
  /* 17   R */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 18   S */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  1,
  /* 19   T */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  1,
  /* 20   U */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,
  /* 21   V */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  1,
  /* 22   W */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,
  /* 23   Y */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,
  /* 24   Z */  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,
  /* 25   * */  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,
  /* 26   X */  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
  };
  // clang-format on


  AminoacidCode toIupac(char aa) {
    const auto it = charToAminoacid.find(aa);
    if (it == charToAminoacid.end()) {
      throw ErrorAminoacidInvalid(aa);
    }
    return it->second;
  }

}// namespace

int lookupAaMatchScoreFromCodes(AminoacidCode x, AminoacidCode y) {
  return scoringMatrixAa[x * AminoacidCode::SIZE + y];//NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
}

int lookupAaMatchScore(char x, char y) {
  return lookupAaMatchScoreFromCodes(toIupac(x), toIupac(y));
}
