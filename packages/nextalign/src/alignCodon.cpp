#include "alignCodon.h"

#include "alignPairwise.h"
#include "helpers.h"
#include "matchNuc.h"


CodonAlignmentResult alignCodon(const std::string& queryPeptide, const std::string& refPeptide) {
  const auto alignment = alignPairwise(queryPeptide, refPeptide, &lookupNucMatchScore, 100);
  return {
    .query = alignment.query,
    .ref = alignment.ref,
    .alignmentScore = alignment.alignmentScore,
  };
}
