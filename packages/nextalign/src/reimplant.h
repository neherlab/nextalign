#pragma once


#include <nextalign/types.h>

#include "alignPairwise.h"

struct AlignmentImproved;
struct CodonAlignmentResult;
struct Gene;


AlignmentImproved reimplant(
  AlignmentImproved& alignmentImproved, const CodonAlignmentResult& codonAlignmentResult, const Gene& gene);
