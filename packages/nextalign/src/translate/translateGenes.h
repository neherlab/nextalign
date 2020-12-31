#pragma once

#include <nextalign/nextalign.h>

#include <vector>

#include "../nextalign_private.h"

struct PeptidesInternal {
  std::vector<PeptideInternal> queryPeptides;
  std::vector<PeptideInternal> refPeptides;
};

PeptidesInternal translateGenes(  //
  const NucleotideSequence& query,//
  const NucleotideSequence& ref,  //
  const GeneMap& geneMap,         //
  const NextalignOptions& options //
);
