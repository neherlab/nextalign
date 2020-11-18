#include <nextalign/types.h>

#include <string>
#include <tuple>

#include "helpers.h"

std::tuple<std::string, GeneMap> parseGb(const std::string& gbContent) {
  NA_UNUSED(gbContent);

  const std::string ref("TODO");
  const GeneMap geneMap;

  return std::make_tuple(ref, geneMap);
}
