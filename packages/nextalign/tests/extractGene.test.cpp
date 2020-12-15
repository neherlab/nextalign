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
    .start = 1,
    .end = 17,
    .strand = "+",
    .frame = 0,
    .length = 16,
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
    .start = 1,
    .end = 17,
    .strand = "+",
    .frame = 0,
    .length = 16,
  };

  const std::string gene_query = extractGeneQuery(query, gene, coordMap);


  EXPECT_EQ(gene_query, expected_gene_query);
}


//TEST(extractGeneQuery, ExtractsQueryGene2) {
//  // i_ref       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24
//  // ref         #  #  #  #  #  #  #  h  e  l  l  o  -  -  i  t  s  -  -  r  e  f  #  #
//  // query       *  *  *  *  *  *  *  h  e  l  l  o  ,  m  i  s  s  i  s  r  e  f  *  *
//  // coordMap    0  1  2  3  4  5  6  7  8  9 10 11       14 15    17 18 19 20 21 22 23 24
//  // i_coordMap  0  1  2  3  4  5  6  7  8  9 10 11       12 13    14 15 16 17 18 19 20 21
//  // gene_ref                         <----------------------------------->
//  // gene_query                       <-------------------->
//
//  const std::string query = "*******hello,missisref**";
//
//  const std::vector<int> coordMap = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24};
//
//  const Gene gene = {
//    .geneName = "Hello",
//    .start = 7,
//    .end = 19,
//    .strand = "+",
//    .frame = 0,
//    .length = 12,
//  };
//
//  const std::string result = extractGeneQuery(query, gene, coordMap);
//  const std::string expected = "heisref";
//
//  EXPECT_EQ(result, expected);
//}
