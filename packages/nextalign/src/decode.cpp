#include "decode.h"

#include <fmt/format.h>
#include <frozen/map.h>
#include <frozen/set.h>
#include <frozen/string.h>
#include <translate.h>

#include "utils/contract.h"


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
  {"TAG", AMINOACID_STOP},
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
  invariant_equal(3, codon.size());

  const auto it = codonTable.find(codon);
  if (it != codonTable.end()) {
    return it->second;
  }

  return 'X';
}


/** A dummy codon. Should never be looked up in the inverse codon table. */
static constexpr const frozen::string DUMMY = "???";

/**
 * Maps each aminoacid to a set of corresponding codons.
 *
 * NOTE: `DUMMY` codons are needed to make set sizes the same, see:
 * https://github.com/serge-sans-paille/frozen/issues/111
 */
static constexpr const frozen::map<char, frozen::set<frozen::string, 6>, 26> inverseCodonTable = {
  {AMINOACID_GAP, {"---", DUMMY, DUMMY, DUMMY, DUMMY, DUMMY}},    //
  {'A', {"GCT", "GCC", "GCA", "GCG", DUMMY, DUMMY}},              //
  {'B', {"AAT", "AAC", "GAT", "GAC", DUMMY, DUMMY}},              //
  {'C', {"TGT", "TGC", DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'D', {"GAT", "GAC", DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'E', {"GAA", "GAG", DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'F', {"TTT", "TTC", DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'G', {"GGT", "GGC", "GGA", "GGG", DUMMY, DUMMY}},              //
  {'H', {"CAT", "CAC", DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'I', {"ATT", "ATC", "ATA", DUMMY, DUMMY, DUMMY}},              //
  /* J  = L | I  */                                               //
  {'K', {"AAA", "AAG", DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'L', {"CTT", "CTC", "CTA", "CTG", "TTA", "TTG"}},              //
  {'M', {"ATG", DUMMY, DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'N', {"AAT", "AAC", DUMMY, DUMMY, DUMMY, DUMMY}},              //
  /* {'O', {"TAG"}}, // Pyrrolysine */                            //
  {'P', {"CCT", "CCC", "CCA", "CCG", DUMMY, DUMMY}},              //
  {'Q', {"CAA", "CAG", DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'R', {"CGT", "CGC", "CGA", "CGG", "AGA", "AGG"}},              //
  {'S', {"TCT", "TCC", "TCA", "TCG", "AGT", "AGC"}},              //
  {'T', {"ACT", "ACC", "ACA", "ACG", DUMMY, DUMMY}},              //
  /* {'U', {"TGA"}}, // Selenocysteine */                         //
  {'V', {"GTT", "GTC", "GTA", "GTG", DUMMY, DUMMY}},              //
  {'W', {"TGG", DUMMY, DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'X', {"---", DUMMY, DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'Y', {"TAT", "TAC", DUMMY, DUMMY, DUMMY, DUMMY}},              //
  {'Z', {"CAA", "CAG", "GAA", "GAG", DUMMY, DUMMY}}, /*=  E | Q *///
  {AMINOACID_START, {"ATG", DUMMY, DUMMY, DUMMY, DUMMY, DUMMY}},  //
  {AMINOACID_STOP, {"TAA", "TGA", "TAG", DUMMY, DUMMY, DUMMY}}    //
};


/**
 * Maps *some* of the aminoacids to a set of corresponding compressed codons.
 *
 * Some of the inverse codons can be compressed using IUPAC nucleotide notation.
 * Not all codons have a compressed form. In this case aminoacid is not present.
 * To be used together with the non-compressed table.
 *
 * NOTE: `DUMMY` codons are needed to make set sizes the same, see:
 * https://github.com/serge-sans-paille/frozen/issues/111
 */
static constexpr const frozen::map<char, frozen::set<frozen::string, 4>, 22> inverseCodonTableCompressed = {
  {'A', {"GCN", DUMMY, DUMMY, DUMMY}},           //
  {'B', {"RAY", DUMMY, DUMMY, DUMMY}},           //
  {'C', {"TGY", DUMMY, DUMMY, DUMMY}},           //
  {'D', {"GAY", DUMMY, DUMMY, DUMMY}},           //
  {'E', {"GAR", DUMMY, DUMMY, DUMMY}},           //
  {'F', {"TTY", DUMMY, DUMMY, DUMMY}},           //
  {'G', {"GGN", DUMMY, DUMMY, DUMMY}},           //
  {'H', {"CAY", DUMMY, DUMMY, DUMMY}},           //
  {'I', {"ATH", DUMMY, DUMMY, DUMMY}},           //
  /* 'J' */                                      //
  {'K', {"AAR", DUMMY, DUMMY, DUMMY}},           //
  {'L', {"CTN", "TTR", "CTY", "YTR"}},           //
  /* 'M' */                                      //
  {'N', {"AAY", DUMMY, DUMMY, DUMMY}},           //
  /* 'O' */                                      //
  {'P', {"CCN", DUMMY, DUMMY, DUMMY}},           //
  {'Q', {"CAR", DUMMY, DUMMY, DUMMY}},           //
  {'R', {"CGN", "AGR", "CGY", "MGR"}},           //
  {'S', {"TCN", "AGY", DUMMY, DUMMY}},           //
  {'T', {"ACN", DUMMY, DUMMY, DUMMY}},           //
  /* 'U' */                                      //
  {'V', {"GTN", DUMMY, DUMMY, DUMMY}},           //
  /* 'W' */                                      //
  {'X', {"---", DUMMY, DUMMY, DUMMY}},           //
  {'Y', {"TAY", DUMMY, DUMMY, DUMMY}},           //
  {'Z', {"SAR", DUMMY, DUMMY, DUMMY}},           //
  /* AMINOACID_START */                          //
  {AMINOACID_STOP, {"TRA", "TAR", DUMMY, DUMMY}},//
};


class ErrorEncodeUnknownAminoacid : public std::runtime_error {
public:
  explicit ErrorEncodeUnknownAminoacid(const Aminoacid& aa)
      : std::runtime_error(fmt::format("When encoding: unknown aminoacid\"{:c}\"", aa)) {}
};


const frozen::set<frozen::string, 6>& getCodons(const Aminoacid& aa) {
  const auto it = inverseCodonTable.find(aa);
  if (it != inverseCodonTable.end()) {
    return it->second;
  }

  throw ErrorEncodeUnknownAminoacid(aa);
}

std::string encode(Aminoacid aa, const std::string_view& hint) {
  precondition_equal(hint.length(), 3);

  const auto& codons = getCodons(aa);

  // Exact match
  const auto& foundExact = codons.find(hint);
  if (foundExact != codons.end()) {
    return std::string{foundExact->data()};
  }

  // Matching 2 first nucleotides
  const auto& found2 = std::find_if(codons.cbegin(), codons.cend(),
    [&hint](const auto& candidate) { return candidate[0] == hint[0] && candidate[1] == hint[1]; });
  if (found2 != codons.end()) {
    return std::string{found2->data()};
  }

  return std::string{DUMMY.data()};
}
