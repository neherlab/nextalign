#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include <nextalign/parseSequences.h>

#include <sstream>

using ExpectedResults = std::vector<std::pair<std::string, std::string>>;

TEST(parseSequences, SanitizesSequences) {
  std::stringstream input;

  // clang-format off
  input << R"(  >  Hello/Sequence/ID1234

 ABC?D
EF.GH

  L*M#N!OP


X Y:)Z
  )"  ;
  // clang-format on

  const auto results = parseSequences(input);

  const ExpectedResults expected = {
    {"Hello/Sequence/ID1234", "ABC?DEF.GHL*MNOPXYZ"},
  };
  EXPECT_EQ(results.size(), expected.size());
  EXPECT_THAT(results, testing::UnorderedElementsAreArray(expected));
}


TEST(parseSequences, ConvertsSequencesToUppercase) {
  std::stringstream input;

  input << R"(>Some/Sequence
  Hi,
  Can you make it uppercase please?
  Cheers!
  )";

  const auto results = parseSequences(input);

  const ExpectedResults expected = {
    {"Some/Sequence", "HICANYOUMAKEITUPPERCASEPLEASE?CHEERS"},
  };
  EXPECT_EQ(results.size(), expected.size());
  EXPECT_THAT(results, testing::UnorderedElementsAreArray(expected));
}

TEST(parseSequences, DeduplicatesSequenceNames) {
  std::stringstream input;

  input << R"(
    >Hello
    ABCD

    >World
    EFGH

    >Foo
    Bar

    >World
    IJKLMN

    >Hello
    OPQRS
  )";

  const auto results = parseSequences(input);

  const ExpectedResults expected = {
    {"Hello", "ABCD"},
    {"World", "EFGH"},
    {"Foo", "BAR"},
    {"World (1)", "IJKLMN"},
    {"Hello (1)", "OPQRS"},
  };
  EXPECT_EQ(results.size(), expected.size());
  EXPECT_THAT(results, testing::UnorderedElementsAreArray(expected));
}

TEST(parseSequences, AssignsSequenceNameToUntitledSequences) {
  std::stringstream input;

  input << R"(
    >
    ABCD
    >
    EFGH
  )";

  const auto results = parseSequences(input);

  const ExpectedResults expected = {
    {"Untitled", "ABCD"},
    {"Untitled (1)", "EFGH"},
  };
  EXPECT_EQ(results.size(), expected.size());
  EXPECT_THAT(results, testing::UnorderedElementsAreArray(expected));
}
