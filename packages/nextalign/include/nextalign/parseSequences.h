#pragma once

#include <istream>
#include <map>
#include <memory>


class FastaStream {
public:
  /** Checks that the stream is in valid state and the next element can be retrieved from it */
  [[nodiscard]] virtual bool good() const = 0;

  /** Retrieves the next sequence in the stream */
  virtual std::pair<std::string, std::string> next() = 0;
};

/** Creates an instance of fasta stream, given a file or string stream */
std::unique_ptr<FastaStream> makeFastaStream(std::istream &istream);

/** Parses all sequences of a given file of string stream */
std::map<std::string, std::string> parseSequences(std::istream &istream);
