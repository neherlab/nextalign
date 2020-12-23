#include "../src/align/alignPairwise.h"

#include <gtest/gtest.h>

#include <string>

#include "../src/match/matchNuc.h"

const int min_length = 5;

TEST(alignPairwise, MatchesIdenticalStrings) {
  std::stringstream input;

  // clang-format off
  const auto qry = toNucleotideSequence("ACGCTCGCT");
  const auto ref = toNucleotideSequence("ACGCTCGCT");
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(ref, result.ref);
  EXPECT_EQ(qry, result.query);
}

TEST(alignPairwise, PadsMissingLeft) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(  "CTCGCT"     );
  const auto ref =    toNucleotideSequence(  "ACGCTCGCT"  );
  const auto qryAln = toNucleotideSequence(  "---CTCGCT"  );
  // FIXME: actual qryAln                    "-C--TCGCT"
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(ref, result.ref);
  EXPECT_EQ(qryAln, result.query);
}

TEST(alignPairwise, PadsMissingLeftSingleMismatch) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(     "GTCCAGCC"  );
  const auto ref =    toNucleotideSequence(  "AGGTACAACCAGCC"  );
  const auto qryAln = toNucleotideSequence(  "--GT----CCAGCC"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(toString(ref), toString(result.ref));
  EXPECT_EQ(toString(qryAln), toString(result.query));
}


TEST(alignPairwise, PadsMissingLeftSingle) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(   "CGCTCGCT"  );
  const auto ref =    toNucleotideSequence(  "ACGCTCGCT"  );
  const auto qryAln = toNucleotideSequence(  "-CGCTCGCT"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(toString(ref), toString(result.ref));
  EXPECT_EQ(toString(qryAln), toString(result.query));
}

TEST(alignPairwise, PadsMissingLeftMismatch) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(      "TGTTACCTGCGC" );
  const auto ref =    toNucleotideSequence( "AAGGTTTATACCTGCGC" );
  const auto qryAln = toNucleotideSequence( "--TGTT---ACCTGCGC" );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(toString(ref), toString(result.ref));
  EXPECT_EQ(toString(qryAln), toString(result.query));
}

TEST(alignPairwise, PadsMissingRight) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(  "ACGCTC"     );
  const auto ref =    toNucleotideSequence(  "ACGCTCGCT"  );
  const auto qryAln = toNucleotideSequence(  "ACGCTC---"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(ref, result.ref);
  EXPECT_EQ(qryAln, result.query);
}

TEST(alignPairwise, HandlesQueryContainedInRef) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(  "ACGCTC"        );
  const auto ref =    toNucleotideSequence(  "GCCACGCTCGCT"  );
  const auto qryAln = toNucleotideSequence(  "---ACGCTC---"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(ref, result.ref);
  EXPECT_EQ(qryAln, result.query);
}

TEST(alignPairwise, HandlesRefContainedInQuery) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(  "GCCACGCTCGCT"  );
  const auto ref =    toNucleotideSequence(     "ACGCTC"     );
  const auto refAln = toNucleotideSequence(  "---ACGCTC---"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(refAln, result.ref);
  EXPECT_EQ(qry, result.query);
}

TEST(alignPairwise, AddsGapsWhenOneMismatch) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(  "GCCACTCCCT"    );
  const auto ref =    toNucleotideSequence(  "GCCACGCTCGCT"  );
  const auto qryAln = toNucleotideSequence(  "GCCAC--TCCCT"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(ref, result.ref);
  EXPECT_EQ(qryAln, result.query);
}

TEST(alignPairwise, AddsGapsInRefWhenOneAmbiguousButMatchingChar) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(  "GCCACGCTCRCT"  );
  const auto ref =    toNucleotideSequence(  "GCCACTCGCT"    );
  const auto refAln = toNucleotideSequence(  "GCCAC--TCGCT"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(refAln, result.ref);
  EXPECT_EQ(qry, result.query);
}

TEST(alignPairwise, CorrectlyAlignsAmbiguousGapPlacingCase) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(  "ACATCTT"        );
  const auto ref =    toNucleotideSequence(  "ACATATGGCACTT"  );
  const auto qryAln = toNucleotideSequence(  "ACAT------CTT"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(ref, result.ref);
  EXPECT_EQ(qryAln, result.query);
}


TEST(alignPairwise, CorrectlyAlignsAmbiguousGapPlacingCaseReversed) {
  std::stringstream input;

  // clang-format off
  const auto qry =    toNucleotideSequence(  "ACATATGGCACTT"  );
  const auto ref =    toNucleotideSequence(  "ACATCTT"        );
  const auto refAln = toNucleotideSequence(  "ACAT------CTT"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(refAln, result.ref);
  EXPECT_EQ(qry, result.query);
}

TEST(alignPairwise, CorrectlyAlignsLongComplexQuery) {
  std::stringstream input;

  // clang-format off
  const auto ref =    toNucleotideSequence(  "CTTGGAGGTTCCGTGGCTAGATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGCGACACTTCGCATGGTGAGCCTTTGT"         );
  const auto qry =    toNucleotideSequence(  "CTTGGAGGTTCCGTGGCTATAAAGATAACAGAACATTCTTGGAATGCTGATCAAGCTCATGGGACANNNNNCATGGTGGACAGCCTTTGT"       );
  const auto refAln = toNucleotideSequence(  "CTTGGAGGTTCCGTGGCTA----GATAACAGAACATTCTTGGAATGCTGATCTTTATAAGCTCATGCGACACTTCGCATGGTG---AGCCTTTGT"  );
  const auto qryAln = toNucleotideSequence(  "CTTGGAGGTTCCGTGGCTATAAAGATAACAGAACATTCTTGGAATGCTGATC-----AAGCTCATGGGACANNNNNCATGGTGGACAGCCTTTGT"  );
  // clang-format on

  const auto result = alignPairwise(qry, ref, min_length);
  EXPECT_EQ(refAln, result.ref);
  EXPECT_EQ(qryAln, result.query);
}
