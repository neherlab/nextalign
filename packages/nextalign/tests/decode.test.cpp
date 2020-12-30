#include "../src/translate/decode.h"

#include <gtest/gtest.h>

#include "../src/alphabet/aminoacids.h"
#include "../src/alphabet/nucleotides.h"

TEST(decode, DecodesGap) {

  EXPECT_EQ(toChar(decode(toNucleotideSequence("---"))), '-');
}

TEST(decode, DecodesValidAminoacid) {
  EXPECT_EQ(toChar(decode(toNucleotideSequence("ATG"))), 'M');
}

TEST(decode, DecodesUnknownToX) {
  EXPECT_EQ(toChar(decode(toNucleotideSequence("HI!"))), 'X');
}
