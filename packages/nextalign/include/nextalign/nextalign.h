#pragma once

#include <istream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>


struct AlgorithmInput {
  int index;
  std::string seqName;
  std::string seq;
};

struct NextalignOptions {
  std::set<std::string> genes;
};

struct Gene {
  std::string geneName;
  int start;
  int end;
  std::string strand;
  int frame;
  int length;
};


using GeneMap = std::map<std::string, Gene>;

struct Alignment {
  std::string query;
  int alignmentScore;
};

struct Insertion {
  int begin;
  int end;
  std::string seq;
};

struct AlignmentImproved : public Alignment {
  std::vector<Insertion> insertions;
};

struct AlgorithmOutput {
  int index;
  std::string seqName;
  bool hasError;
  AlignmentImproved result;
  std::exception_ptr error;
};


AlignmentImproved nextalign(
  const std::string& query, const std::string& ref, const GeneMap& geneMap, const NextalignOptions& options);

/**
 * Parses genemap in GFF format from a file or string stream
 *
 * @see GFF format reference at https://www.ensembl.org/info/website/upload/gff.html
 */
GeneMap parseGeneMapGff(std::istream& is);

class FastaStream {
public:
  FastaStream() = default;

  virtual ~FastaStream() = default;

  FastaStream(const FastaStream& other) = delete;

  FastaStream& operator=(const FastaStream& other) = delete;

  FastaStream(FastaStream&& other) = delete;

  FastaStream& operator=(const FastaStream&& other) = delete;

  /** Checks that the stream is in valid state and the next element can be retrieved from it */
  [[nodiscard]] virtual bool good() const = 0;

  /** Retrieves the next sequence in the stream */
  virtual AlgorithmInput next() = 0;
};

/** Creates an instance of fasta stream, given a file or string stream */
std::unique_ptr<FastaStream> makeFastaStream(std::istream& istream);

/** Parses all sequences of a given file or string stream */
std::vector<AlgorithmInput> parseSequences(std::istream& istream);
