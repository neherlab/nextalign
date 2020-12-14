#pragma once

#include <string_view>
#include <vector>

struct Gene;

std::string_view extractGeneRef(const std::string_view& ref, const Gene& gene);

std::string_view extractGeneQuery(const std::string_view& query, const Gene& gene, const std::vector<int>& coordMap);
