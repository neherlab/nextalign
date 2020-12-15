#pragma once

#include <fmt/format.h>

#include <exception>
#include <string_view>
#include <vector>

struct Gene;

class ErrorExtractGeneLengthInvalid : public std::runtime_error {
public:
  explicit ErrorExtractGeneLengthInvalid(const std::string& gene, int numGaps)
      : std::runtime_error(fmt::format("When extracting a gene: genes expected to have a number of deletions that is a "
                                       "multiple of 3, but Gene \"{:s}\" has {:d}",
          gene, numGaps)) {}
};

std::string_view extractGeneRef(const std::string_view& ref, const Gene& gene);

std::string extractGeneQuery(const std::string_view& query, const Gene& gene, const std::vector<int>& coordMap);
