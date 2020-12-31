#include <nextalign/nextalign.h>

#include <string>

#include "align/alignPairwise.h"
#include "align/getGapOpenCloseScores.h"
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

  // TODO: hoist this out of the loop
  const auto gapOpenCloseNuc = getGapOpenCloseScoresCodonAware(ref, geneMap, options);
  const auto gapOpenCloseAA = getGapOpenCloseScoresFlat(ref, options);

  const auto alignment = alignPairwise(query, ref, gapOpenCloseNuc, 100);

  std::vector<Peptide> queryPeptides;
  std::vector<Peptide> refPeptides;
  if (!options.genes.empty()) {
    const auto peptidesInternal = translateGenes(alignment.query, alignment.ref, geneMap, gapOpenCloseAA, options);
    queryPeptides = map(peptidesInternal.queryPeptides, std::function<Peptide(PeptideInternal)>(toPeptideExternal));
    refPeptides = map(peptidesInternal.refPeptides, std::function<Peptide(PeptideInternal)>(toPeptideExternal));
  }

  const auto stripped = stripInsertions(alignment.ref, alignment.query);

  NextalignResult result;
  result.query = toString(stripped.queryStripped);
  result.alignmentScore = alignment.alignmentScore;
  result.refPeptides = refPeptides;
  result.queryPeptides = queryPeptides;
  result.insertions = map(stripped.insertions, std::function<Insertion(InsertionInternal)>(toInsertionExternal));

  return result;
}
