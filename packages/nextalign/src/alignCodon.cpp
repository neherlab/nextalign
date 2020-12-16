#include "alignCodon.h"

#include "alignPairwise.h"
#include "helpers.h"
#include "matchAa.h"


CodonAlignmentResult alignCodon(const std::string& queryPeptide, const std::string& refPeptide) {
  const auto alignment = alignPairwise(queryPeptide, refPeptide, &lookupAaMatchScore, 10);

  //  postcondition_equal(queryPeptide.size(), alignment.query.size());
  //  postcondition_equal(refPeptide.size(), alignment.ref.size());

  return {
    .queryPeptide = alignment.query,
    .refPeptide = alignment.ref,
    .alignmentScore = alignment.alignmentScore,
  };
}
