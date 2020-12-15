#pragma once

#include <nextalign/types.h>

#include <string>


struct NextalignOptions;

AlignmentImproved nextalign(
  const std::string& query, const std::string& ref, const GeneMap& geneMap, const NextalignOptions& options);
