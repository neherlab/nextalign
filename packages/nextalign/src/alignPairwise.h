#pragma once

#include <nextalign/types.h>

#include <string>

struct NextalignOptions;

Alignment alignPairwise(const std::string& query, const std::string& ref, const NextalignOptions& options);
