#include <nextalign/nextalign.h>

#include <string>

#include "align/alignPairwise.h"
#include "alphabet/nucleotides.h"
#include "strip/stripInsertions.h"
#include "utils/contract.h"
#include "utils/map.h"
#include "utils/safe_cast.h"


Insertion toExternal(const InsertionInternal& ins) {
  return Insertion{.begin = ins.begin, .end = ins.end, .seq = toString(ins.seq)};
}

AlignmentImproved nextalign(const NucleotideSequence& query, const NucleotideSequence& ref, const GeneMap& geneMap,
  const NextalignOptions& options) {

  const auto alignment = alignPairwise(query, ref, 100);
  const auto stripped = stripInsertions(alignment.ref, alignment.query);

  AlignmentImproved result;
  result.query = toString(stripped.queryStripped);
  result.alignmentScore = alignment.alignmentScore;
  result.insertions = map(stripped.insertions, toExternal, std::vector<Insertion>{});

  return result;
}
