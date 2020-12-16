#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <array>

#include "../src/mapCoordinates.h"

TEST(mapCoordinates, MapsSimple) {
  const std::string ref = "ACTC---CGTG---A";
  constexpr const auto expected = std::array<int, 9>{0, 1, 2, 3, 7, 8, 9, 10, 14};
  EXPECT_THAT(mapCoordinates(ref), testing::ElementsAreArray(expected));
}
