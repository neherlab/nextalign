#include "alignCodon.h"

#include "alignPairwise.h"
#include "helpers.h"
#include "matchAa.h"


CodonAlignmentResult alignCodon(const std::string& queryPeptide, const std::string& refPeptide) {
  const auto alignment = alignPairwise(queryPeptide, refPeptide, &lookupAaMatchScore, 100);
  return {
    .queryPeptide = alignment.query,
    .refPeptide = alignment.ref,
    .alignmentScore = alignment.alignmentScore,
  };
}
