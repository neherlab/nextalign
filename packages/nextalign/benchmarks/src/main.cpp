#include <benchmark/benchmark.h>

#include "utils/getData.h"

const auto [sequences, ref, geneMap, totalNucs] = getData();

// clang-format off
#include "Nextalign.benchmark.h"
#include "AlignPairwise.benchmark.h"
#include "SeedMatching.benchmark.h"
// clang-format on


BENCHMARK_MAIN();
