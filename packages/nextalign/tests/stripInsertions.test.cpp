#include <gtest/gtest.h>

#include "../src/strip/stripInsertions.h"

TEST(stripInsertions, StripsAnInsertion) {
  // clang-format off
  const auto reference = toNucleotideSequence("ACT---CTCTACTCTAC");
  const auto query     = toNucleotideSequence("ACTGCGCTCTAC---AC");
  const auto expected  = toNucleotideSequence("ACTCTCTAC---AC");
  // clang-format on

  const auto& res = stripInsertions(reference, query);

  EXPECT_EQ(res.queryStripped, expected);
  EXPECT_EQ(res.insertions[0].begin, 3);
  EXPECT_EQ(res.insertions[0].end, 6);
  EXPECT_EQ(res.insertions[0].seq, toNucleotideSequence("GCG"));
}
