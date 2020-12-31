#pragma once

#include <nextalign/nextalign.h>

#include <gsl/string_span>

template<typename Letter>
using SequenceSpan = gsl::basic_string_span<Letter, gsl::dynamic_extent>;


struct PeptideInternal {
  std::string name;
  AminoacidSequence seq;
};
