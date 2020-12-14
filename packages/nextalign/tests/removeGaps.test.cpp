
#include "../src/removeGaps.h"

#include <gtest/gtest.h>

constexpr const auto* const INPUT = "--MY-SPACEBAR---IS--BROKEN-SEND---HELP-";
constexpr const auto* const OUTPUT = "MYSPACEBARISBROKENSENDHELP";


TEST(removeGaps, RemovesGaps) {
  const std::string input(INPUT);
  EXPECT_EQ(removeGaps(input), OUTPUT);
  EXPECT_EQ(input, INPUT);
}

TEST(removeGapsInPlace, RemovesGapsInPlace) {
  std::string input(INPUT);
  removeGapsInPlace(input);
  EXPECT_EQ(input, OUTPUT);
}
