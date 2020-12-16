#include "stripInsertions.h"

#include <nextalign/types.h>

#include <string>

#include "safeCast.h"
#include "utils/contract.h"


StripInsertionsResult stripInsertions(const std::string& ref, const std::string& query) {
  const int refLength = safe_cast<int>(ref.size());
  const int queryLength = safe_cast<int>(query.size());
  precondition_equal(refLength, queryLength);

  StripInsertionsResult result;
  result.queryStripped.reserve(refLength);

  int insertionStart = -1;
  std::string currentInsertion;
  for (int i = 0; i < refLength; ++i) {
    const auto& c = ref[i];
    if (c == '-') {
      if (currentInsertion.empty()) {
        currentInsertion = query[i];
        insertionStart = i;
      } else {
        currentInsertion += query[i];
      }
    } else {
      result.queryStripped += query[i];
      if (!currentInsertion.empty()) {
        const auto length = safe_cast<int>(currentInsertion.size());
        const auto end = insertionStart + length;

        result.insertions.emplace_back(Insertion{.begin = insertionStart, .end = end, .seq = currentInsertion});

        currentInsertion = "";
        insertionStart = -1;
      }
    }
  }

  precondition_less_equal(result.queryStripped.size(), refLength);
  return result;
}
