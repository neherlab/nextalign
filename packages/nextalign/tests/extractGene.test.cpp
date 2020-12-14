#include "../src/extractGene.h"

#include <gtest/gtest.h>
#include <nextalign/types.h>

#include <string>
#include <string_view>
#include <vector>

#include "../src/mapCoordinates.h"


TEST(extractGeneRef, ExtractsRefGene) {
  /*
    gene:                |               |
    ref:                ACTCCG---TGCGCTGCTG
    query:              ACTGGTCTAGTCG---TAG
    extractGene(ref) =  'CTCCGTGCGCTGC'
    extractGene(query) ='CTGGTCTAGTCGT'
   */

  const std::string ref = "ACTCCG---TGCGCTGCTG";

  const Gene gene = {
    .geneName = "Hello",
    .start = 1,
    .end = 17,
    .strand = "+",
    .frame = 0,
    .length = 16,
  };

  const std::string result = extractGeneRef(ref, gene);
  const std::string expected = "CTCCGTGCGCTGC";

  EXPECT_EQ(result, expected);
}


TEST(extractGeneQuery, ExtractsQueryGene) {
  /*
    gene:                |               |
    ref:                ACTCCG---TGCGCTGCTG
    query:              ACTGGTCTAGTCG---TAG
    extractGene(ref) =  'CTCCGTGCGCTGC'
    extractGene(query) ='CTGGTCTAGTCGT'
   */

  const std::string ref = "ACTCCG---TGCGCTGCTG";
  const std::string query = "ACTGGTCTAGTCG---TAG";

  const auto coordMap = mapCoordinates(ref);

  const Gene gene = {
    .geneName = "Hello",
    .start = 1,
    .end = 17,
    .strand = "+",
    .frame = 0,
    .length = 16,
  };

  const std::string result = extractGeneQuery(query, gene, coordMap);
  const std::string expected = "CTGGTCTAGTCGT";

  EXPECT_EQ(result, expected);
}


//TEST(extractGeneQuery, ExtractsQueryGene2) {
//  // i_ref       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
//  // ref         #  #  #  #  #  #  #  h  e  l  l  o  -  -  i  t  -  s  -  -  -  r  e  f  #  #
//  // query       *  *  *  *  *  *  *  h  e  -  -  i  -  s  r  e  f  *  *
//  // coordMap    0  1  2  3  4  5  6  7  8  9 10 11       14 15    17 18 19 20 21 22 23 24 25
//  // i_coordMap  0  1  2  3  4  5  6  7  8  9 10 11       12 13    14 15 16 17 18 19 20 21 22
//  // gene_ref                         <----------------------------------->
//  // gene_query                       <-------------------->
//
//  const std::string query = "*******helloitsref**";
//
//  const std::vector<int> coordMap = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25};
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
//  const std::string_view result = extractGeneQuery(query, gene, coordMap);
//  const std::string expected = "helloits";
//
//  EXPECT_EQ(result, expected);
//}
