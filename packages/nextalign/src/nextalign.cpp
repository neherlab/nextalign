

// * parse and prepare reference and genemap
//* stream the sequences
//    - for each sequence
//       * seed matching
//       * pairwise align (banded)
//
//       * new stuff:
//         * for each gene
//          - extract sequence
//          - translate
//          - codon align
//          - reverse translate
//          - reimplant
//       * output

#include <fmt/format.h>

#include <string>

#include "alignPairwise.h"
#include "helpers.h"
#include "nextalign/types.h"
#include "parseGb.h"

class ErrorGeneMapGeneNotFound : std::runtime_error {
  std::string formatError(const std::string& geneName) const {
    return fmt::format("Error: gene '{:s}' not found in gene map", geneName);
  }

public:
  ErrorGeneMapGeneNotFound(const std::string& geneName) : std::runtime_error(formatError(geneName)) {}
};

void matchSeeds() {}

template<typename Container, typename Element>
bool includes(const Container& container, const Element& element) {
  return container.find(element) != container.end();
}


const std::string extractGeneSequence(const std::string& ref, const Gene& gene) {
  return "TODO";
}

Alignment improveAlignment(const Alignment& alignment, const GeneMap& geneMap, const NextalignOptions& options) {

  for (const auto& geneName : options.genes) {

    const auto& found = geneMap.find(geneName);
    if (found != geneMap.end()) {
      throw ErrorGeneMapGeneNotFound(geneName);
    }
  }

  // For each gene of the virus
  for (const auto& [geneName, gene] : geneMap) {
    // Ignore gene, if it's not included in the subset of genes to be used
    if (!includes(options.genes, geneName)) {
      continue;
    }

    //
    const auto& geneSeq = extractGeneSequence(alignment.ref, gene);
  }

  // * for each gene
  //  - extract sequence
  //  - translate
  //  - codon align
  //  - reverse translate
  //  - reimplant

  const Alignment alignmentImproved = alignment;

  return alignmentImproved;
}

Alignment nextalign(
  const std::string& query, const std::string& ref, const GeneMap& geneMap, const NextalignOptions& options) {
  matchSeeds();

  const auto minimalLength = 0;
  const auto alignment = alignPairwise(query, ref, options);

  const auto alignmentImproved = improveAlignment(alignment, geneMap, options);
  return alignmentImproved;
}
