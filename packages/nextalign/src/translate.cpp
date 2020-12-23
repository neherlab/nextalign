#include "translate.h"

#include <string>
#include <string_view>

#include "decode.h"
#include "safeCast.h"
#include "utils/contract.h"


Peptide translate(const std::string_view& seq) {
  precondition_divisible_by(seq.size(), 3);

  const int seqLength = safe_cast<int>(seq.size());
  const int peptideLength = seqLength / 3;


  Peptide peptide(peptideLength, AMINOACID_GAP);
  for (int i_aa = 0; i_aa < peptideLength; ++i_aa) {
    const auto i_nuc = i_aa * 3;
    const auto codon = seq.substr(i_nuc, 3);
    const auto aminoacid = decode(codon);

    invariant_less(i_aa, peptide.size());
    peptide[i_aa] = aminoacid;
  }

  postcondition_equal(peptide.size(), peptideLength);
  postcondition_equal(seq.size(), peptide.size() * 3);
  return peptide;
}