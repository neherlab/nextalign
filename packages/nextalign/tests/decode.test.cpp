#include <gtest/gtest.h>

#include "../src/decode.h"

TEST(decode, DecodesGap) {
  EXPECT_EQ(decode("---"), AMINOACID_GAP);
}

TEST(decode, DecodesValidAminoacid) {
  EXPECT_EQ(decode("ATG"), 'M');
}

TEST(decode, DecodesUnknownToX) {
  EXPECT_EQ(decode("HI!"), AMINOACID_UNKNOWN);
}
