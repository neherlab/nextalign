#include "reimplant.h"

#include <nextalign/types.h>

#include <string>
#include <string_view>

#include "alignCodon.h"
#include "alignPairwise.h"
#include "decode.h"
#include "helpers.h"
#include "safeCast.h"
#include "utils/contract.h"


Alignment reimplant(Alignment& alignmentImproved, const CodonAlignmentResult& codonAlignmentResult, const Gene& gene) {
  auto& query = alignmentImproved.query;
  const auto queryLength = safe_cast<int>(query.size());

  const auto& refPeptide = std::string_view{codonAlignmentResult.refPeptide};

  const auto& queryPeptide = std::string_view{codonAlignmentResult.queryPeptide};
  const auto queryPeptideLength = safe_cast<int>(queryPeptide.size());

  for (int i_aa = 0; i_aa < queryPeptideLength; ++i_aa) {
    const int i_nuc = gene.start + i_aa * 3;

    invariant_less(i_aa, refPeptide.size());
    invariant_less(i_aa, queryPeptide.size());
    invariant_less(i_nuc, query.size());

    const auto& aa_query = queryPeptide[i_aa];
    if (aa_query == AMINOACID_GAP) {
      query[i_nuc + 0] = '-';
      query[i_nuc + 1] = '-';
      query[i_nuc + 2] = '-';
    }
  }

  postcondition_equal(query.size(), queryLength);// query length should not change
  return alignmentImproved;
}
