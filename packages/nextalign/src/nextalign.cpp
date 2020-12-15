#include <fmt/format.h>
#include <nextalign/types.h>
#include <utils/contract.h>

#include <string>

#include "alignCodon.h"
#include "alignPairwise.h"
#include "extractGene.h"
#include "helpers.h"
#include "mapCoordinates.h"
#include "reimplant.h"
#include "safeCast.h"
#include "translate.h"

class ErrorGeneMapGeneNotFound : std::runtime_error {
  static std::string formatError(const std::string& geneName) {
    return fmt::format("Error: gene '{:s}' not found in gene map", geneName);
  }

public:
  explicit ErrorGeneMapGeneNotFound(const std::string& geneName) : std::runtime_error(formatError(geneName)) {}
};

void matchSeeds() {}


Alignment alignBetter(const Alignment& alignment, const GeneMap& geneMap, const NextalignOptions& options) {
  Alignment alignmentImproved = alignment;

  const auto& ref = alignment.ref;
  const auto& query = alignment.query;
  //  const auto coordMap = mapCoordinates(ref);
  const auto coordMap = std::vector<int>{};

  // For each gene in the requested subset
  for (const auto& geneName : options.genes) {
    // TODO: Should probably validate gene names before even running
    const auto& found = geneMap.find(geneName);
    if (found == geneMap.end()) {
      throw ErrorGeneMapGeneNotFound(geneName);
    }

    const auto& gene = found->second;

    const auto& refGene = extractGeneRef(ref, gene);
    const auto refPeptide = translate(refGene);

    const auto& queryGene = extractGeneQuery(query, gene, coordMap);
    const auto queryPeptide = translate(queryGene);

    const CodonAlignmentResult codonAlignmentResult = alignCodon(refPeptide, queryPeptide);

    reimplant(alignmentImproved, codonAlignmentResult, gene);
  }

  return alignmentImproved;
}

Alignment nextalign(
  const std::string& query, const std::string& ref, const GeneMap& geneMap, const NextalignOptions& options) {
  matchSeeds();// TODO: mmmm..?
  const auto alignment = alignPairwise(query, ref, &lookupNucMatchScore, 100);
  auto alignmentImproved = alignBetter(alignment, geneMap, options);
  return alignmentImproved;
}
