#include "decode.h"

#include <frozen/map.h>
#include <frozen/string.h>
#include <translate.h>


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
