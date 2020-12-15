#pragma once

#include <string>

struct CodonAlignmentResult {
  std::string queryPeptide;
  std::string refPeptide;
  int alignmentScore;
};

CodonAlignmentResult alignCodon(const std::string& queryPeptide, const std::string& refPeptide);
