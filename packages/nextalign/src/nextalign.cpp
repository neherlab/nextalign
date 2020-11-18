#include <fmt/format.h>
#include <nextalign/types.h>

#include <string>

#include "alignPairwise.h"
#include "helpers.h"

class ErrorGeneMapGeneNotFound : std::runtime_error {
  static std::string formatError(const std::string& geneName) {
    return fmt::format("Error: gene '{:s}' not found in gene map", geneName);
  }

public:
  explicit ErrorGeneMapGeneNotFound(const std::string& geneName) : std::runtime_error(formatError(geneName)) {}
};

void matchSeeds() {}

std::string extractGeneSequence(const std::string& ref, const Gene& gene) {
  NA_UNUSED(ref);
  NA_UNUSED(gene);
  return "GENE SEQUENCE";
}

std::string translate(const std::string& seq) {
  NA_UNUSED(seq);
  return "PEPTIDE";
}

struct CodonAlignmentResult {
  // SOMETHING
};


CodonAlignmentResult codonAlign(const std::string& peptide) {
  NA_UNUSED(peptide);
  return {};
}

Alignment reimplant(const Alignment& alignment, CodonAlignmentResult codonAlignmentResult) {
  NA_UNUSED(alignment.ref);
  NA_UNUSED(alignment.query);
  NA_UNUSED(alignment.alignmentScore);
  NA_UNUSED(codonAlignmentResult);
  return alignment;
}


Alignment alignBetter(const Alignment& alignment, const GeneMap& geneMap, const NextalignOptions& options) {
  Alignment alignmentImproved = alignment;

  // For each gene in the requested subset
  for (const auto& geneName : options.genes) {
    // TODO: Should probably validate gene names before even running
    const auto& found = geneMap.find(geneName);
    if (found != geneMap.end()) {
      throw ErrorGeneMapGeneNotFound(geneName);
    }

    const auto& gene = found->second;

    const std::string geneSeq = extractGeneSequence(alignment.ref, gene);// TODO: ref or query or both?
    const std::string peptide = translate(geneSeq);

    const CodonAlignmentResult codonAlignmentResult = codonAlign(peptide);

    alignmentImproved = reimplant(alignmentImproved, codonAlignmentResult);
  }

  return alignmentImproved;
}

Alignment nextalign(
  const std::string& query, const std::string& ref, const GeneMap& geneMap, const NextalignOptions& options) {
  matchSeeds();// TODO: mmmm..?
  const auto alignment = alignPairwise(query, ref, options);
  auto alignmentImproved = alignBetter(alignment, geneMap, options);
  return alignmentImproved;
}
