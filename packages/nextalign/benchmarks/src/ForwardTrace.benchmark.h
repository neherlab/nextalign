#pragma once

#include <benchmark/benchmark.h>

#include <numeric>
#include <vector>

#include "../include/nextalign/nextalign.h"
#include "../src/alignPairwise.h"
#include "utils/getData.h"
#include "utils/setCounters.h"


class ForwardTraceBench : public benchmark::Fixture {
protected:
  std::vector<SeedAlignment> seedAlignments;


  ForwardTraceBench() {
    const auto n = NUM_SEQUENCES_AVG;
    seedAlignments.resize(n);
    for (int i = 0; i < n; ++i) {
      const auto& [seqName, query] = sequences[i];
      seedAlignments[i] = seedAlignment(query, reference);
    }
  }
};


/**
 * Average benchmark for scoreMatrix().
 * Runs `scoreMatrix()` for NUM_SEQUENCES_AVG sequences and averages the result.
 * This is an estimate of runtime performance in a real world scenario, when many sequences are ran in a batch.
 */

BENCHMARK_DEFINE_F(ForwardTraceBench, Average)(benchmark::State& st) {
  const auto n = NUM_SEQUENCES_AVG;
  ForwardTrace forwardTrace;
  st.SetComplexityN(totalNucs);

  for (const auto& _ : st) {
    for (int i = 0; i < n; ++i) {
      const auto& [seqName, query] = sequences[i];
      const auto& seedAlignment = seedAlignments[i];
      benchmark::DoNotOptimize(
        forwardTrace = scoreMatrix(query, reference, seedAlignment.bandWidth, seedAlignment.meanShift));
    }
  }

  setCounters(st, n);
}


BENCHMARK_REGISTER_F(ForwardTraceBench, Average)
  ->Unit(benchmark::kMillisecond)//
  ->Complexity(benchmark::oNSquared)
  ->Iterations(3);


/**
 * Variation benchmark for scoreMatrix().
 * Runs `scoreMatrix()` for NUM_SEQUENCES_VAR sequences and shows results per sequence.
 * This shows variation or runtime between different sequences.
 */
BENCHMARK_DEFINE_F(ForwardTraceBench, Variation)(benchmark::State& st) {
  const auto& index = st.range(0);
  const auto& [seqName, query] = sequences[index];
  SeedAlignment seed = seedAlignments[index];
  ForwardTrace forwardTrace;
  st.SetLabel(seqName);
  st.SetComplexityN(query.size());

  for (const auto& _ : st) {
    benchmark::DoNotOptimize(forwardTrace = scoreMatrix(query, reference, seed.bandWidth, seed.meanShift));
  }

  setCounters(st, 1);
}

BENCHMARK_REGISTER_F(ForwardTraceBench, Variation)
  ->DenseRange(0, NUM_SEQUENCES_VAR, 1)//
  ->Unit(benchmark::kMillisecond)      //
  ->Complexity(benchmark::oNSquared);
