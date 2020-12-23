#include <alphabet/nucleotides.h>
#include <fmt/format.h>
#include <nextalign/types.h>
#include <utils/contract.h>

#include <gsl/string_span>
#include <string>

#include "alignPairwise.h"
#include "src/utils/safe_cast.h"
#include "stripInsertions.h"
#include "utils/map.h"

Insertion toExternal(const InsertionInternal& ins) {
  return Insertion{.begin = ins.begin, .end = ins.end, .seq = toString(ins.seq)};
}

AlignmentImproved nextalign(
  const std::string& queryStr, const std::string& refStr, const GeneMap& geneMap, const NextalignOptions& options) {

  const auto query = toNucleotideSequence(queryStr);
  const auto ref = toNucleotideSequence(refStr);

  const auto alignment = alignPairwise(query, ref, 100);
  const auto stripped = stripInsertions(alignment.ref, alignment.query);

  AlignmentImproved result;
  result.query = toString(stripped.queryStripped);
  result.alignmentScore = alignment.alignmentScore;
  result.insertions = map(stripped.insertions, toExternal, std::vector<Insertion>{});

  return result;
}
