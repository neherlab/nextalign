#include <gtest/gtest.h>

#include <string>

#include "../src/align/alignPairwise.h"
#include "../src/align/getGapOpenCloseScores.h"
#include "../src/match/matchNuc.h"

const int min_length = 5;

TEST(alignPairwise, AlignsSimpleOneGene) {
  const NextalignOptions options = {
    .gapOpenInFrame = 1,
    .gapOpenOutOfFrame = 2,
    .genes = {"Gene 1"},
  };

  GeneMap geneMap = {//
    {"Gene 1",       //
      Gene{
        .geneName = "Gene1",
        .start = 2,
        .end = 19,
        .strand = "+",
        .frame = 0,
        .length = 17,
      }}};


  // clang-format off
  // gapOpenCosts                             22221221221221221   221222
  const auto qry =    toNucleotideSequence(  "GCATGAAGAATATCATT" "GCTTTG"  );
  const auto ref =    toNucleotideSequence(  "GCATGAAGAATATCATT" "GCTTTG"  );
  const auto refAln = toNucleotideSequence(  "GCATGAAGAATATCATT---GCTTTG"  );
  const auto qryAln = toNucleotideSequence(  "-CATG---AATATCATTAATGCTTTG"  );
  // clang-format on

  const std::vector<int> gapOpenCosts = getGapOpenCloseScoresCodonAware(ref, geneMap, options);

  const auto result = alignPairwise(qry, ref, gapOpenCosts, min_length);
  EXPECT_EQ(toString(ref), toString(result.ref));
  EXPECT_EQ(toString(qryAln), toString(result.query));
}
