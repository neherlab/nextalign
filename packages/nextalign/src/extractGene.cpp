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


std::string extractGeneRef(const std::string_view& ref, const Gene& gene) {
  precondition_less(gene.length, ref.size());
  precondition_less_equal(gene.length, ref.size());

  const auto unstripped = ref.substr(gene.start, gene.length);
  auto stripped = removeGaps(unstripped);

  // HACK: adjust gene length to be a multiple of 3
  resizeToMultipleOf3InPlace(stripped);

  // Length of the gene should not exceed the length of the sequence
  invariant_less_equal(stripped.size(), ref.size());

  // Gene length should be a multiple of 3
  invariant_divisible_by(stripped.size(), 3);

  return stripped;
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


  const auto unstripped = query.substr(start, length);
  auto stripped = removeGaps(unstripped);

  // HACK: adjust gene length to be a multiple of 3
  resizeToMultipleOf3InPlace(stripped);

  // Length of the gene should not exceed the length of the sequence
  invariant_less_equal(stripped.size(), query.size());

  // Gene length should be a multiple of 3
  invariant_divisible_by(stripped.size(), 3);

  return stripped;
}
