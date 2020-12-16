#pragma once

#include <nextalign/types.h>

#include <gsl/string_span>
#include <string_view>


struct CodonAlignmentResult;

void translateReverseInPlace(                              //
  const std::string_view& queryGene,                       //
  const CodonAlignmentResult& codonAlignmentResult,        //
  /* out */ gsl::string_span<>& reverseTranslatedQueryGene,//
  /* out */ std::vector<Insertion>& insertions             //

);
