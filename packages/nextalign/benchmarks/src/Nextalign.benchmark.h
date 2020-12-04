#pragma once

#include <benchmark/benchmark.h>

#include <numeric>
#include <vector>

#include "../include/nextalign/nextalign.h"
#include "../src/alignPairwise.h"
#include "utils/getData.h"
#include "utils/setCounters.h"


/**
 * Average benchmark for nextalign().
 * Runs `nextalign()` for NUM_SEQUENCES_AVG sequences and averages the result.
 * This is an estimate of runtime performance in a real world scenario, when many sequences are ran in a batch.
 */
void NextalignAverage(benchmark::State& st) {
  const auto n = NUM_SEQUENCES_AVG;
  const NextalignOptions options = {};
  Alignment aln;
  st.SetComplexityN(totalNucs);

  for (const auto& _ : st) {
    for (int i = 0; i < n; ++i) {
      const auto& [seqName, query] = sequences[i];
      benchmark::DoNotOptimize(aln = nextalign(query, reference, geneMap, options));
    }
  }

  setCounters(st, n);
}

BENCHMARK(NextalignAverage)      //
  ->Unit(benchmark::kMillisecond)//
  ->Complexity(benchmark::oNSquared)
  ->Iterations(3);


/**
 * Variation benchmark for nextalign().
 * Runs `nextalign()` for NUM_SEQUENCES_VAR sequences and shows results per sequence.
 * This shows variation or runtime between different sequences.
 */
void NextalignVariation(benchmark::State& st) {
  const auto& index = st.range(0);
  const auto& [seqName, query] = sequences[index];
  const NextalignOptions options = {};
  Alignment aln;
  st.SetLabel(seqName);
  st.SetComplexityN(query.size());

  for (const auto& _ : st) {
    benchmark::DoNotOptimize(aln = nextalign(query, reference, geneMap, options));
  }

  setCounters(st, 1);
}

BENCHMARK(NextalignVariation)          //
  ->DenseRange(0, NUM_SEQUENCES_VAR, 1)//
  ->Unit(benchmark::kMillisecond)      //
  ->Complexity(benchmark::oNSquared);
