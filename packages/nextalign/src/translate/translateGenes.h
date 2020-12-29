#pragma once

#include <nextalign/nextalign.h>

#include <vector>

#include "../nextalign_private.h"

std::vector<PeptideInternal> translateGenes(//
  const NucleotideSequence& query,          //
  const NucleotideSequence& ref,            //
  const GeneMap& geneMap,                   //
  const NextalignOptions& options           //
);
