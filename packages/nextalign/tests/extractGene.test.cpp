#include "../src/extractGene.h"

#include <gtest/gtest.h>
#include <nextalign/types.h>

#include <string>
#include <string_view>
#include <vector>

#include "../src/mapCoordinates.h"


TEST(extractGeneRef, ExtractsRefGene) {
  // clang-format off
  //                                       2             17
  // gene                                  |              |
  const std::string ref =               "ACTCCG---TGCGCTGCTG";
  // query                               ACTGGTCTAGTCG---TAG
  const std::string expected_gene_ref =   "TCCGTGCGCTGC";
  // expected_gene_query                   TGGTCTAGTCGT
  // clang-format on

  const Gene gene = {
    .geneName = "Hello",
    .start = 2,
    .end = 17,
    .strand = "+",
    .frame = 0,
    .length = 15,
  };

  const std::string gene_ref = extractGeneRef(ref, gene);

  EXPECT_EQ(gene_ref, expected_gene_ref);
}


TEST(extractGeneQuery, ExtractsQueryGene) {
  // clang-format off
  //                                       2             17
  // gene                                  |              |
  const std::string ref =               "ACTCCG---TGCGCTGCTG";
  const std::string query =             "ACTGGTCTAGTCG---TAG";
  // expected gene_ref                     TCCGTGCGCTGC
  const std::string expected_gene_query = "TGGTCTAGTCGT";
  // clang-format on

  //  const auto coordMap = mapCoordinates(ref);
  const auto coordMap = std::vector<int>{};

  const Gene gene = {
    .geneName = "Hello",
    .start = 2,
    .end = 17,
    .strand = "+",
    .frame = 0,
    .length = 15,
  };

  const std::string gene_query = extractGeneQuery(query, gene, coordMap);


  EXPECT_EQ(gene_query, expected_gene_query);
}
