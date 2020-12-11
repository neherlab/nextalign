#include "decode.h"

#include <frozen/map.h>
#include <frozen/string.h>
#include <translate.h>



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


static constexpr const frozen::map<frozen::string, char, 65> codonTable = {
  {"---", AMINOACID_GAP},
  {"AAA", 'K'},
  {"AAC", 'N'},
  {"AAG", 'K'},
  {"AAT", 'N'},
  {"ACA", 'T'},
  {"ACC", 'T'},
  {"ACG", 'T'},
  {"ACT", 'T'},
  {"AGA", 'R'},
  {"AGC", 'S'},
  {"AGG", 'R'},
  {"AGT", 'S'},
  {"ATA", 'I'},
  {"ATC", 'I'},
  {"ATG", 'M'},
  {"ATT", 'I'},
  {"CAA", 'Q'},
  {"CAC", 'H'},
  {"CAG", 'Q'},
  {"CAT", 'H'},
  {"CCA", 'P'},
  {"CCC", 'P'},
  {"CCG", 'P'},
  {"CCT", 'P'},
  {"CGA", 'R'},
  {"CGC", 'R'},
  {"CGG", 'R'},
  {"CGT", 'R'},
  {"CTA", 'L'},
  {"CTC", 'L'},
  {"CTG", 'L'},
  {"CTT", 'L'},
  {"GAA", 'E'},
  {"GAC", 'D'},
  {"GAG", 'E'},
  {"GAT", 'D'},
  {"GCA", 'A'},
  {"GCC", 'A'},
  {"GCG", 'A'},
  {"GCT", 'A'},
  {"GGA", 'G'},
  {"GGC", 'G'},
  {"GGG", 'G'},
  {"GGT", 'G'},
  {"GTA", 'V'},
  {"GTC", 'V'},
  {"GTG", 'V'},
  {"GTT", 'V'},
  {"TAA", '*'},
  {"TAC", 'Y'},
  {"TAG", '*'},
  {"TAT", 'Y'},
  {"TCA", 'S'},
  {"TCC", 'S'},
  {"TCG", 'S'},
  {"TCT", 'S'},
  {"TGA", '*'},
  {"TGC", 'C'},
  {"TGG", 'W'},
  {"TGT", 'C'},
  {"TTA", 'L'},
  {"TTC", 'F'},
  {"TTG", 'L'},
  {"TTT", 'F'},
};

Aminoacid decode(const std::string_view& codon) {
  const auto it = codonTable.find(codon);
  if (it != codonTable.end()) {
    return it->second;
  }

  return 'X';
}
