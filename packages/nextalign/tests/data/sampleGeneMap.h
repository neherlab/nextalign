#pragma once

#include <nextalign/types.h>

#include <vector>

const std::vector<Gene> sampleGeneArray = {//
  Gene{
    .geneName = "E",
    .start = 26245,
    .end = 26472,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "M",
    .start = 26523,
    .end = 27191,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "N",
    .start = 28274,
    .end = 29533,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF10",
    .start = 29558,
    .end = 29674,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF14",
    .start = 28734,
    .end = 28955,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF1a",
    .start = 266,
    .end = 13468,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF1b",
    .start = 13468,
    .end = 21555,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF3a",
    .start = 25393,
    .end = 26220,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF6",
    .start = 27202,
    .end = 27387,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF7a",
    .start = 27394,
    .end = 27759,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF7b",
    .start = 27756,
    .end = 27887,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF8",
    .start = 27894,
    .end = 28259,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "ORF9b",
    .start = 28284,
    .end = 28577,
    .strand = "+",
    .frame = 0,
  },
  Gene{
    .geneName = "S",
    .start = 21563,
    .end = 25384,
    .strand = "+",
    .frame = 0,
  }};


GeneMap makeGeneMap(const std::vector<Gene>& genes) {
  GeneMap geneMap;
  std::transform(genes.begin(), genes.end(), std::inserter(geneMap, geneMap.end()),
    [](const Gene& gene) { return std::make_pair(gene.geneName, gene); });
  return geneMap;
}

const GeneMap sampleGeneMap = makeGeneMap(sampleGeneArray);
