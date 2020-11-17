#include "alignPairwise.h"

#include <string>


Alignment alignPairwise(const std::string& query, const std::string& ref, int minimalLength) {
  (void) minimalLength;// unused

  const auto alignedQuery = query;
  const auto alignedRef = ref;
  const auto alignmentScore = 0;

  return {.query = alignedQuery, .ref = alignedRef, .alignmentScore = alignmentScore};
}
