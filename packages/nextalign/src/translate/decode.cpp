#include "decode.h"

#include <fmt/format.h>
#include <frozen/map.h>
#include <frozen/set.h>
#include <frozen/string.h>

#include "utils/contract.h"


//using NucleotideSequenceFrozen = frozen::basic_string<Nucleotide>;

//#define make_codon(c1, c2, c3) NucleotideSequenceFrozen({Nucleotide::c1, Nucleotide::c2, Nucleotide::c3})

//static constexpr const frozen::map<NucleotideSequenceFrozen, Aminoacid, 65> codonTable = {
//  {make_codon(GAP, GAP, GAP), Aminoacid::GAP},
//  {make_codon(A, A, A), Aminoacid::K},
//  {make_codon(A, A, C), Aminoacid::N},
//  {make_codon(A, A, G), Aminoacid::K},
//  {make_codon(A, A, T), Aminoacid::N},
//  {make_codon(A, C, A), Aminoacid::T},
//  {make_codon(A, C, C), Aminoacid::T},
//  {make_codon(A, C, G), Aminoacid::T},
//  {make_codon(A, C, T), Aminoacid::T},
//  {make_codon(A, G, A), Aminoacid::R},
//  {make_codon(A, G, C), Aminoacid::S},
//  {make_codon(A, G, G), Aminoacid::R},
//  {make_codon(A, G, T), Aminoacid::S},
//  {make_codon(A, T, A), Aminoacid::I},
//  {make_codon(A, T, C), Aminoacid::I},
//  {make_codon(A, T, G), Aminoacid::M},
//  {make_codon(A, T, T), Aminoacid::I},
//  {make_codon(C, A, A), Aminoacid::Q},
//  {make_codon(C, A, C), Aminoacid::H},
//  {make_codon(C, A, G), Aminoacid::Q},
//  {make_codon(C, A, T), Aminoacid::H},
//  {make_codon(C, C, A), Aminoacid::P},
//  {make_codon(C, C, C), Aminoacid::P},
//  {make_codon(C, C, G), Aminoacid::P},
//  {make_codon(C, C, T), Aminoacid::P},
//  {make_codon(C, G, A), Aminoacid::R},
//  {make_codon(C, G, C), Aminoacid::R},
//  {make_codon(C, G, G), Aminoacid::R},
//  {make_codon(C, G, T), Aminoacid::R},
//  {make_codon(C, T, A), Aminoacid::L},
//  {make_codon(C, T, C), Aminoacid::L},
//  {make_codon(C, T, G), Aminoacid::L},
//  {make_codon(C, T, T), Aminoacid::L},
//  {make_codon(G, A, A), Aminoacid::E},
//  {make_codon(G, A, C), Aminoacid::D},
//  {make_codon(G, A, G), Aminoacid::E},
//  {make_codon(G, A, T), Aminoacid::D},
//  {make_codon(G, C, A), Aminoacid::A},
//  {make_codon(G, C, C), Aminoacid::A},
//  {make_codon(G, C, G), Aminoacid::A},
//  {make_codon(G, C, T), Aminoacid::A},
//  {make_codon(G, G, A), Aminoacid::G},
//  {make_codon(G, G, C), Aminoacid::G},
//  {make_codon(G, G, G), Aminoacid::G},
//  {make_codon(G, G, T), Aminoacid::G},
//  {make_codon(G, T, A), Aminoacid::V},
//  {make_codon(G, T, C), Aminoacid::V},
//  {make_codon(G, T, G), Aminoacid::V},
//  {make_codon(G, T, T), Aminoacid::V},
//  {make_codon(T, A, A), Aminoacid::STOP},
//  {make_codon(T, A, C), Aminoacid::Y},
//  {make_codon(T, A, G), Aminoacid::STOP},
//  {make_codon(T, A, T), Aminoacid::Y},
//  {make_codon(T, C, A), Aminoacid::S},
//  {make_codon(T, C, C), Aminoacid::S},
//  {make_codon(T, C, G), Aminoacid::S},
//  {make_codon(T, C, T), Aminoacid::S},
//  {make_codon(T, G, A), Aminoacid::STOP},
//  {make_codon(T, G, C), Aminoacid::C},
//  {make_codon(T, G, G), Aminoacid::W},
//  {make_codon(T, G, T), Aminoacid::C},
//  {make_codon(T, T, A), Aminoacid::L},
//  {make_codon(T, T, C), Aminoacid::F},
//  {make_codon(T, T, G), Aminoacid::L},
//  {make_codon(T, T, T), Aminoacid::F},
//};

Aminoacid decode(const NucleotideSequenceView& codon) {
//  invariant_equal(3, codon.size());
//
//  const auto it = codonTable.find(codon);
//  if (it != codonTable.end()) {
//    return it->second;
//  }

  return Aminoacid::X;
}
