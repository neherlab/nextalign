#include <string>

#include "alignPairwise.h"
#include "helpers.h"
#include "nextalign/types.h"
#include "seqan2.h"

Alignment alignPairwise(const std::string& query, const std::string& ref, const NextalignOptions& options) {
  NA_UNUSED(options);


  const AlignOptionsOverlap alignOptions = {
    .band = 5,
    .score_match = 1,
    .score_mismatch = -1,
    .score_gapext = 1,
    .score_gapopen = -1,
    .cut_flanks = 0,
  };

  const AlignmentResult& result = align_overlap(query, ref, alignOptions);

  return {.query = result.seq1, .ref = result.seq2, .alignmentScore = result.score};
}
