#include <gtest/gtest.h>

#include "../src/stripInsertions.h"

TEST(stripInsertions, StripsAnInsertion) {
  // clang-format off
  const std::string reference = "ACT---CTCTACTCTAC";
  const std::string query     = "ACTGCGCTCTAC---AC";
  const std::string expected  = "ACTCTCTAC---AC";
  // clang-format on

  const auto& res = stripInsertions(reference, query);

  EXPECT_EQ(res.queryStripped, expected);
  EXPECT_EQ(res.insertions[0].begin, 3);
  EXPECT_EQ(res.insertions[0].end, 6);
  EXPECT_EQ(res.insertions[0].seq, "GCG");
}
