#pragma once

#include <benchmark/benchmark.h>

#include <numeric>
#include <vector>

#include "../include/nextalign/nextalign.h"
#include "../src/alignPairwise.h"
#include "utils/getData.h"
#include "utils/setCounters.h"


class BackwardTraceBench : public benchmark::Fixture {
protected:
  std::vector<SeedAlignment> seedAlignments;
  std::vector<ForwardTrace> forwardTraces;

  BackwardTraceBench() {
    const auto n = NUM_SEQUENCES_AVG;
    seedAlignments.resize(n);
    forwardTraces.resize(n);
    for (int i = 0; i < n; ++i) {
      const auto& [seqName, query] = sequences[i];
      seedAlignments[i] = seedAlignment(query, ref);
      forwardTraces[i] = scoreMatrix(query, ref, seedAlignments[i].bandWidth, seedAlignments[i].meanShift);
    }
  }
};


/**
 * Average benchmark for backTrace().
 * Runs `backTrace()` for NUM_SEQUENCES_AVG sequences and averages the result.
 * This is an estimate of runtime performance in a real world scenario, when many sequences are ran in a batch.
 */

BENCHMARK_DEFINE_F(BackwardTraceBench, Average)(benchmark::State& st) {
  const auto n = NUM_SEQUENCES_AVG;
  Alignment aln;
  st.SetComplexityN(totalNucs);

  for (const auto& _ : st) {
    for (int i = 0; i < n; ++i) {
      const auto& [seqName, query] = sequences[i];
      const auto& seed = seedAlignments[i];
      const auto& forwardTrace = forwardTraces[i];
      benchmark::DoNotOptimize(aln = backTrace(query, ref, forwardTrace.scores, forwardTrace.paths, seed.meanShift));
    }
  }

  setCounters(st, n);
}


BENCHMARK_REGISTER_F(BackwardTraceBench, Average)
  ->Unit(benchmark::kMillisecond)//
  ->Complexity(benchmark::oNSquared)
  ->Iterations(3);


class BackwardTraceBench2 : public benchmark::Fixture {
protected:
  std::vector<SeedAlignment> seedAlignments;
  std::vector<ForwardTrace> forwardTraces;

  BackwardTraceBench2() {
    const auto n = NUM_SEQUENCES_AVG;
    seedAlignments.resize(n);
    forwardTraces.resize(n);
    for (int i = 0; i < n; ++i) {
      const auto& [seqName, query] = sequences[i];
      seedAlignments[i] = seedAlignment(query, ref);
      forwardTraces[i] = scoreMatrix(query, ref, seedAlignments[i].bandWidth, seedAlignments[i].meanShift);
    }
  }
};

/**
 * Variation benchmark for backTrace().
 * Runs `backTrace()` for NUM_SEQUENCES_VAR sequences and shows results per sequence.
 * This shows variation or runtime between different sequences.
 */
BENCHMARK_DEFINE_F(BackwardTraceBench2, Variation)(benchmark::State& st) {
  const auto& index = st.range(0);
  const auto& [seqName, query] = sequences[index];
  SeedAlignment seed = seedAlignments[index];
  ForwardTrace forwardTrace = forwardTraces[index];
  Alignment aln;
  st.SetLabel(seqName);
  st.SetComplexityN(query.size());

  for (const auto& _ : st) {
    benchmark::DoNotOptimize(aln = backTrace(query, ref, forwardTrace.scores, forwardTrace.paths, seed.meanShift));
  }

  setCounters(st, 1);
}

BENCHMARK_REGISTER_F(BackwardTraceBench2, Variation)
  ->DenseRange(0, NUM_SEQUENCES_VAR, 1)//
  ->Unit(benchmark::kMillisecond)      //
  ->Complexity(benchmark::oNSquared)
  ->Iterations(0);
