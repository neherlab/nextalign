#include <gtest/gtest.h>

#include <string>

#include "../src/align/alignPairwise.h"
#include "../src/align/getGapOpenCloseScores.h"
#include "../src/match/matchNuc.h"

const int min_length = 5;

TEST(alignPairwise, AlignsSimpleOneGene) {
  const NextalignOptions options = {
    .gapOpenInFrame = -5,
    .gapOpenOutOfFrame = -6,
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
  const auto ref =    toNucleotideSequence(  "GCATGAGGAATCTCAGTGCTTTG"  );
  const auto refAln = toNucleotideSequence(  "GCATGAGGAATCTCAGT---GCTTTG"  );
  const auto qryAln = toNucleotideSequence(  "-CATG---AATCTCAGTAATGCTTTG"  );
  const auto qry =    toNucleotideSequence(   "CATGAATCTCAGTAATGCTTTG"  );
  // clang-format on

  const std::vector<int> gapOpenCosts = getGapOpenCloseScoresCodonAware(ref, geneMap, options);

  const auto result = alignPairwise(qry, ref, gapOpenCosts, min_length);
  EXPECT_EQ(19*3-5-5, result.alignmentScore);
  EXPECT_EQ(toString(refAln), toString(result.ref));
  EXPECT_EQ(toString(qryAln), toString(result.query));
}
