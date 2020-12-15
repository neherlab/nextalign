#include "extractGene.h"

#include <nextalign/types.h>

#include <string_view>
#include <vector>

#include "alignPairwise.h"
#include "removeGaps.h"
#include "safeCast.h"

void resizeToMultipleOf3InPlace(std::string& seq) {
  const auto length = safe_cast<int>(seq.size());
  const auto excess = length % 3;
  seq.resize(length - excess);
}


std::string_view extractGeneRef(const std::string_view& ref, const Gene& gene) {
  precondition_less(gene.length, ref.size());
  precondition_less_equal(gene.length, ref.size());

  // TODO: should we also strip gaps from the raw gene?
  return ref.substr(gene.start, gene.length);
}


/**
 * Extracts gene from the query sequence according to coordinate map relative to the reference sequence
 */
std::string extractGeneQuery(const std::string_view& query, const Gene& gene, const std::vector<int>& coordMap) {
  precondition_less(gene.start, coordMap.size());
  precondition_less(gene.end, coordMap.size());

  // Transform gene coordinates according to coordinate map
  const auto start = coordMap[gene.start];
  const auto end = coordMap[gene.end];// TODO: `gene.end` -1 or not?
  const auto length = end - start;

  // Start and end should be within bounds
  invariant_less(start, query.size());
  invariant_less(end, query.size());

  const auto unstripped = query.substr(start, length);
  auto stripped = removeGaps(unstripped);

  const auto numGaps = safe_cast<int>(unstripped.size() - stripped.size());
  if (numGaps % 3 != 0) {
    throw ErrorExtractGeneLengthInvalid(gene.geneName, numGaps);
  }

  invariant_less_equal(stripped.size(), query.size());// Length of the gene should not exceed the length of the sequence
  invariant_divisible_by(stripped.size(), 3);         // Gene length should be a multiple of 3
  return stripped;
}
