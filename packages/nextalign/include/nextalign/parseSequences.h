#pragma once

#include <nextalign/types.h>

#include <istream>
#include <map>
#include <memory>

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
