#pragma once

#include <string>

struct CodonAlignmentResult {
  std::string query;
  std::string ref;
  int alignmentScore;
};

CodonAlignmentResult alignCodon(const std::string& queryPeptide, const std::string& refPeptide);
