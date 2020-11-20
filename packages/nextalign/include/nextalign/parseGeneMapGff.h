#pragma once

#include <nextalign/types.h>

#include <sstream>

/**
 * Parses genemap in GFF format from a file or string stream
 *
 * @see GFF format reference at https://www.ensembl.org/info/website/upload/gff.html
 */
GeneMap parseGeneMapGff(std::istream& is);
