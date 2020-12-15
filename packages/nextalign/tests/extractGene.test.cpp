#include "../src/extractGene.h"

#include <gtest/gtest.h>
#include <nextalign/types.h>

#include <string>
#include <string_view>
#include <vector>

#include "../src/mapCoordinates.h"


TEST(extractGeneRef, ExtractsRefGene) {
  // clang-format off
  // gene                                    2          14
  // gene                                    |           |
  const std::string ref =                 "ACTCCGTGCGCTGCTG";
  const std::string ref_aligned =         "ACTCCG---TGCGCTGCTG";
  const std::string query_aligned =       "ACTGGTCTA---GCA";
  const std::string expected_gene_ref =     "TCCGTGCGCTGC";
  const std::string expected_gene_query =   "TGGTCTAGC";
  // clang-format on

  const Gene gene = {
    .geneName = "Hello",
    .start = 2,
    .end = 14,
    .strand = "+",
    .frame = 0,
    .length = 12,
  };

  const auto gene_ref = std::string(extractGeneRef(ref, gene));

  EXPECT_EQ(gene_ref, expected_gene_ref);
}


TEST(extractGeneQuery, ExtractsQueryGene) {
  // clang-format off
  // gene                                    2          14
  // gene                                    |           |
  const std::string ref =                 "ACTCCGTGCGCTGCTG";
  const std::string ref_aligned =         "ACTCCG---TGCGCTGCTG";
  const std::string query_aligned =       "ACTGGTCTA---GCA";
  const std::string expected_gene_ref =     "TCCGTGCGCTGC";
  const std::string expected_gene_query =   "TGGTCTAGC";
  // clang-format on

  const Gene gene = {
    .geneName = "Hello",
    .start = 2,
    .end = 14,
    .strand = "+",
    .frame = 0,
    .length = 12,
  };

  const auto coordMap = mapCoordinates(ref);
  const std::string gene_query = extractGeneQuery(query_aligned, gene, coordMap);

  EXPECT_EQ(gene_query, expected_gene_query);
}
