#include "extractGene.h"

#include <nextalign/types.h>

#include <string_view>
#include <vector>

#include "alignPairwise.h"
#include "removeGaps.h"

std::string extractGeneRef(const std::string_view& ref, const Gene& gene) {
  precondition_less(gene.length, ref.size());
  precondition_less_equal(gene.length, ref.size());
  return removeGaps(ref.substr(gene.start, gene.length));
}

/**
 * Extracts gene from the query sequence according to coordinate map relative to the reference sequence
 */
std::string extractGeneQuery(const std::string_view& query, const Gene& gene, const std::vector<int>& coordMap) {
  //  precondition_less(gene.start, coordMap.size());
  //  precondition_less(gene.end, coordMap.size());

  const auto start = gene.start;
  const auto end = gene.end;
  const auto length = end - start;

  //  // Transform gene coordinates according to coordinate map
  //  const auto start = coordMap[gene.start];
  //  const auto end = coordMap[gene.end];// TODO: `gene.end` -1 or not?
  //  const auto length = end - start;

  // Start and end should be within bounds
  invariant_less(start, query.size());
  invariant_less(end, query.size());

  // Length of the gene should not exceed the length of the sequence
  invariant_less_equal(length, query.size());

  return removeGaps(query.substr(start, length));
}
