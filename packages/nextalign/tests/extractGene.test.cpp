#include "../src/extractGene.h"

#include <gtest/gtest.h>
#include <nextalign/types.h>

#include <string>
#include <string_view>
#include <vector>

TEST(extractGeneQuery, ExtractsQueryGene) {
  // i_ref       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
  // ref         #  #  #  #  #  #  #  h  e  l  l  o  -  -  i  t  -  s  -  -  -  r  e  f  #  #
  // query       *  *  *  *  *  *  *  h  e  l  l  o  i  t  s  r  e  f  *  *
  // coordMap    0  1  2  3  4  5  6  7  8  9 10 11       14 15    17 18 19 20 21 22 23 24 25
  // i_coordMap  0  1  2  3  4  5  6  7  8  9 10 11       12 13    14 15 16 17 18 19 20 21 22
  // gene_ref                         <----------------------------------->
  // gene_query                       <-------------------->

  const std::string query = "*******helloitsref**";

  const std::vector<int> coordMap = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 17, 18, 19, 20, 21, 22, 23, 24, 25};

  const Gene gene = {
    .geneName = "Hello",
    .start = 7,
    .end = 19,
    .strand = "+",
    .frame = 0,
    .length = 12,
  };

  const std::string_view result = extractGeneQuery(query, gene, coordMap);
  const std::string expected = "helloits";

  EXPECT_EQ(result, expected);
}
