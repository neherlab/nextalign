#include "translateReverse.h"

#include <nextalign/types.h>

#include <string>
#include <string_view>

#include "alignCodon.h"
#include "alignPairwise.h"
#include "decode.h"
#include "helpers.h"
#include "safeCast.h"
#include "utils/contract.h"


void translateReverseInPlace(                              //
  const std::string_view& queryGene,                       //
  const CodonAlignmentResult& codonAlignmentResult,        //
  /* out */ gsl::string_span<>& reverseTranslatedQueryGene,//
  /* out */ std::vector<Insertion>& insertions             //

) {
  auto& query = reverseTranslatedQueryGene;
  const auto queryLength = safe_cast<int>(query.size());

  const auto& refPeptide = std::string_view{codonAlignmentResult.refPeptide};

  const auto& queryPeptide = std::string_view{codonAlignmentResult.queryPeptide};
  const auto queryPeptideLength = safe_cast<int>(queryPeptide.size());

  for (int i_aa = 0; i_aa < queryPeptideLength; ++i_aa) {
    const int i_nuc = i_aa * 3;

    invariant_less(i_aa, refPeptide.size());
    invariant_less(i_aa, queryPeptide.size());
    invariant_less(i_nuc, query.size());

    const auto& aa_query = queryPeptide[i_aa];
    if (aa_query == AMINOACID_GAP) {
      query[i_nuc + 0] = '-';
      query[i_nuc + 1] = '-';
      query[i_nuc + 2] = '-';
    } else {
      const auto& seq = queryGene.substr(i_nuc, 3);
      query[i_nuc + 0] = seq[0];
      query[i_nuc + 1] = seq[1];
      query[i_nuc + 2] = seq[2];
    }

    const auto& aa_ref = refPeptide[i_aa];
    invariant(!((aa_query == AMINOACID_GAP) && (aa_ref == AMINOACID_GAP)));
    if (aa_ref == AMINOACID_GAP) {
      const auto& seq = queryGene.substr(i_nuc, 3);
      insertions.emplace_back(Insertion{.begin = i_nuc, .end = i_nuc + 3, .seq = std::string{seq}});
    }
  }

  postcondition_equal(query.size(), queryLength);// query length should not change
}
