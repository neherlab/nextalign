#include "../src/alignPairwise.h"

#include <benchmark/benchmark.h>

#include <fstream>
#include <map>
#include <numeric>
#include <vector>

#include "../include/nextalign/nextalign.h"
#include "../include/nextalign/parseGeneMapGff.h"
#include "../include/nextalign/parseSequences.h"

constexpr const int NUM_SEQUENCES_VAR = 10;// Number of sequences to process per "Variation" benchmark
constexpr const int NUM_SEQUENCES_AVG = 30;// Number of sequences to process per "Average" benchmark

inline void setCounters(benchmark::State& st, int nSeq) {
  // How many sequences can be processed per second
  st.counters["seq/s"] = benchmark::Counter(st.iterations() * nSeq, benchmark::Counter::kIsRate);

  // How many milliseconds it takes to process one sequence
  st.counters["ms/seq"] =
    benchmark::Counter(st.iterations() * nSeq / 1000.0, benchmark::Counter::kIsRate | benchmark::Counter::kInvert);
}

auto getData() {
  std::vector<std::pair<std::string, std::string>> sequences;
  std::ifstream fastaFile("../../../../../data/example/sequences.fasta");
  const auto sequencesMap = parseSequences(fastaFile);
  std::copy(sequencesMap.cbegin(), sequencesMap.cend(), back_inserter(sequences));

  assert(sequences.size() < NUM_SEQUENCES_VAR);
  const auto totalNucs =
    std::accumulate(sequences.cbegin(), sequences.cbegin() + NUM_SEQUENCES_VAR, 0, [](int total, const auto& seq) {
      total += seq.second.size();
      return total;
    });

  std::ifstream refFile("../../../../../data/example/reference.txt");
  const auto refSeqs = parseSequences(refFile);
  const auto ref = refSeqs.begin()->second;

  std::ifstream genemapFile("../../../../../data/example/genemap.gff");
  const auto geneMap = parseGeneMapGff(genemapFile);

  return std::make_tuple(sequences, ref, geneMap, totalNucs);
}

const auto [sequences, ref, geneMap, totalNucs] = getData();

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
      benchmark::DoNotOptimize(aln = nextalign(query, ref, geneMap, options));
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
 * Runs `nextalign()` for NUM_SEQUENCES sequences and shows results per sequence.
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
    benchmark::DoNotOptimize(aln = nextalign(query, ref, geneMap, options));
  }

  setCounters(st, 1);
}

BENCHMARK(NextalignVariation)          //
  ->DenseRange(0, NUM_SEQUENCES_VAR, 1)//
  ->Unit(benchmark::kMillisecond)      //
  ->Complexity(benchmark::oNSquared)
  ->Iterations(0);


void AlignPairwiseAverage(benchmark::State& st) {
  const auto n = NUM_SEQUENCES_AVG;
  const NextalignOptions options = {};
  Alignment aln;
  st.SetComplexityN(totalNucs);

  for (const auto& _ : st) {
    for (int i = 0; i < n; ++i) {
      const auto& [seqName, query] = sequences[i];
      benchmark::DoNotOptimize(aln = alignPairwise(query, ref, 100));
    }
  }

  setCounters(st, n);
}

BENCHMARK(AlignPairwiseAverage)  //
  ->Unit(benchmark::kMillisecond)//
  ->Complexity(benchmark::oNSquared)
  ->Iterations(3);


void AlignPairwiseVariation(benchmark::State& st) {
  const auto& index = st.range(0);
  const auto& [seqName, query] = sequences[index];
  Alignment aln;
  st.SetLabel(seqName);
  st.SetComplexityN(query.size());

  for (const auto& _ : st) {
    benchmark::DoNotOptimize(aln = alignPairwise(query, ref, 100));
  }

  setCounters(st, 1);
}

BENCHMARK(AlignPairwiseVariation)      //
  ->DenseRange(0, NUM_SEQUENCES_VAR, 1)//
  ->Unit(benchmark::kMillisecond)      //
  ->Complexity(benchmark::oNSquared)
  ->Iterations(0);


BENCHMARK_MAIN();
