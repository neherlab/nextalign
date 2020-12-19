#include <fmt/format.h>
#include <nextalign/types.h>
#include <utils/contract.h>

#include <gsl/string_span>
#include <string>

#include "alignCodon.h"
#include "alignPairwise.h"
#include "extractGene.h"
#include "helpers.h"
#include "mapCoordinates.h"
#include "safeCast.h"
#include "stripInsertions.h"
#include "translate.h"
#include "translateReverse.h"

class ErrorGeneMapGeneNotFound : std::runtime_error {
  static std::string formatError(const std::string& geneName) {
    return fmt::format("Error: gene '{:s}' not found in gene map", geneName);
  }

public:
  explicit ErrorGeneMapGeneNotFound(const std::string& geneName) : std::runtime_error(formatError(geneName)) {}
};


AlignmentImproved alignBetter(
  const std::string& ref, const Alignment& alignment, const GeneMap& geneMap, const NextalignOptions& options) {

  std::string newQueryMemory(ref.size(), '-');
  gsl::string_span<> newQuery{newQueryMemory};

  std::string newRefMemory(ref.size(), '-');
  gsl::string_span<> newRef{newRefMemory};

  AlignmentImproved alignmentImproved = {
    // base Alignment
    {
      .query = newQueryMemory,
      .ref = newRefMemory,
      .alignmentScore = alignment.alignmentScore,
    },
    {}// insertions
  };

  const auto coordMap = mapCoordinates(alignment.ref);

  // Each position in the raw ref sequence should have a corresponding mapped position in aligned ref sequence
  invariant_equal(coordMap.size(), ref.size());

  std::vector<Insertion> insertions;

  // For each gene in the requested subset
  for (const auto& geneName : options.genes) {
    // TODO: Should probably validate gene names before even running
    const auto& found = geneMap.find(geneName);
    if (found == geneMap.end()) {
      throw ErrorGeneMapGeneNotFound(geneName);
    }

    const auto& gene = found->second;

    // TODO: can be done once during initialization
    const auto& refGene = extractGeneRef(ref, gene);
    const auto refPeptide = translate(refGene);

    const auto& queryGene = extractGeneQuery(alignment.query, gene, coordMap);
    const auto queryPeptide = translate(queryGene);

    const CodonAlignmentResult codonAlignmentResult = alignCodon(refPeptide, queryPeptide);

    // TODO: Can it be that `gene.length` is different from `length` here?
    //  If not, this can be removed and we can use `gene.length`.
    //  But make sure to not overflow/underflow the buffer.
    const int start = coordMap[gene.start];
    const int end = coordMap[gene.end];
    const int length = end - start;
    auto reverseTranslatedQueryGene = newQuery.subspan(start, length);
    translateReverseInPlace(
      queryGene, codonAlignmentResult, /* out */ reverseTranslatedQueryGene, /* out */ insertions);

    //    auto reverseTranslatedRefGene = newRef.subspan(gene.start, gene.length);
    //    std::vector<Insertion> ignored;
    //    translateReverseInPlace(refGene, codonAlignmentResult, /* out */ reverseTranslatedRefGene, /* out */ ignored);

    alignmentImproved.alignmentScore += codonAlignmentResult.alignmentScore;
  }

  return alignmentImproved;
}

//TODO: extract constants
std::vector<int> gapOpenCloseScores(
  const std::string& ref, const GeneMap& geneMap, const NextalignOptions& options) {
  std::vector<int> gapOpenClose(ref.size());
  fill(gapOpenClose.begin(), gapOpenClose.end(), -2);

  for (const auto& geneName : options.genes) {
    // TODO: Should probably validate gene names before even running
    const auto& found = geneMap.find(geneName);
    if (found == geneMap.end()) {
      throw ErrorGeneMapGeneNotFound(geneName);
    }

    const auto& gene = found->second;
    for (int i=gene.start; i<=gene.end; i+=3){
      gapOpenClose[i] = -1;
    }
  }
  return gapOpenClose;
}

AlignmentImproved nextalign(
  const std::string& query, const std::string& ref, const GeneMap& geneMap, const NextalignOptions& options) {

  // TODO: This wants to be done once upstream (doesn't depend on query)
  std::vector<int> gapOpenClose = gapOpenCloseScores(ref, geneMap, options);
  const auto alignment = alignPairwise(query, ref, &lookupNucMatchScore, gapOpenClose, 100);

  //  const auto alignmentBetter = alignBetter(ref, alignment, geneMap, options);

  const auto stripped = stripInsertions(alignment.ref, alignment.query);

  AlignmentImproved result;
  result.query = stripped.queryStripped;
  result.ref = alignment.ref;
  result.alignmentScore = alignment.alignmentScore;
  result.insertions = stripped.insertions;

  return result;
}
