#pragma once

#include <nextalign/types.h>

#include <string>
#include <tuple>


std::tuple<std::string, GeneMap> parseGb(const std::string& gbContent);
