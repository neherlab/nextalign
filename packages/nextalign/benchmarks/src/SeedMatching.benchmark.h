#pragma once

#include <benchmark/benchmark.h>

#include <numeric>
#include <vector>

#include "../include/nextalign/nextalign.h"
#include "../src/alignPairwise.cpp"
#include "utils/getData.h"
#include "utils/setCounters.h"


/**
 * Average benchmark for seedAlignment().
 * Runs `seedAlignment()` for NUM_SEQUENCES_AVG sequences and averages the result.
 * This is an estimate of runtime performance in a real world scenario, when many sequences are ran in a batch.
 */
void SeedMatchAverage(benchmark::State& st) {
  const auto n = NUM_SEQUENCES_AVG;
  SeedAlignment aln;
  st.SetComplexityN(totalNucs);

  for (const auto& _ : st) {
    for (int i = 0; i < n; ++i) {
      const auto& [seqName, query] = sequences[i];
      benchmark::DoNotOptimize(aln = seedAlignment(query, ref));
    }
  }

  setCounters(st, n);
}

BENCHMARK(SeedMatchAverage)      //
  ->Unit(benchmark::kMillisecond)//
  ->Complexity(benchmark::oNSquared)
  ->Iterations(3);


/**
 * Variation benchmark for seedAlignment().
 * Runs `seedAlignment()` for NUM_SEQUENCES_VAR sequences and shows results per sequence.
 * This shows variation or runtime between different sequences.
 */
void SeedMatchVariation(benchmark::State& st) {
  const auto& index = st.range(0);
  const auto& [seqName, query] = sequences[index];
  SeedAlignment aln;
  st.SetLabel(seqName);
  st.SetComplexityN(query.size());

  for (const auto& _ : st) {
    benchmark::DoNotOptimize(aln = seedAlignment(query, ref));
  }

  setCounters(st, 1);
}

BENCHMARK(SeedMatchVariation)          //
  ->DenseRange(0, NUM_SEQUENCES_VAR, 1)//
  ->Unit(benchmark::kMillisecond)      //
  ->Complexity(benchmark::oNSquared)
  ->Iterations(0);
