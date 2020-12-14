#include "mapCoordinates.h"

#include <string_view>
#include <vector>

#include "safeCast.h"

/**
 * Makes a map from reference coordinates to alignment coordinates
 * (Excluding the gaps in reference).
 *
 * Example:
 *   012345678901234
 *   ACTC---CGTG---A -> aln_coord = [0,1,2,3,7,8,9,10,14]
 */
std::vector<int> mapCoordinates(const std::string_view& ref) {
  const auto refLength = safe_cast<int>(ref.size());

  std::vector<int> coordMap;
  coordMap.reserve(refLength);
  for (int i = 0; i < refLength; ++i) {
    if (ref[i] != '-') {
      coordMap.push_back(i);
    }
  }

  coordMap.shrink_to_fit();
  return coordMap;
}
