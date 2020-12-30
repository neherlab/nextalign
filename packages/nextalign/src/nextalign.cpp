#include <nextalign/nextalign.h>

#include <string>

#include "align/alignPairwise.h"
#include "alphabet/nucleotides.h"
#include "strip/stripInsertions.h"
#include "translate/translateGenes.h"
#include "utils/contract.h"
#include "utils/map.h"
#include "utils/safe_cast.h"

Peptide toPeptideExternal(const PeptideInternal& peptide) {
  return Peptide{.name = peptide.name, .seq = toString(peptide.seq)};
}

Insertion toInsertionExternal(const InsertionInternal& ins) {
  return Insertion{.begin = ins.begin, .end = ins.end, .seq = toString(ins.seq)};
}

NextalignResult nextalign(const NucleotideSequence& query, const NucleotideSequence& ref, const GeneMap& geneMap,
  const NextalignOptions& options) {

  const auto alignment = alignPairwise(query, ref, 100);

  const auto queryPeptides = translateGenes(alignment.query, ref, geneMap, options);


  const auto stripped = stripInsertions(alignment.ref, alignment.query);

  NextalignResult result;
  result.query = toString(stripped.queryStripped);
  result.alignmentScore = alignment.alignmentScore;
  result.queryPeptides = map(queryPeptides, std::function<Peptide(PeptideInternal)>(toPeptideExternal));
  result.insertions = map(stripped.insertions, std::function<Insertion(InsertionInternal)>(toInsertionExternal));

  return result;
}
