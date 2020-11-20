#include <fmt/format.h>
#include <nextalign/parseGeneMapGff.h>

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/range/algorithm_ext/erase.hpp>
#include <istream>
#include <string>

#include "csvParser.h"


constexpr const auto GFF_NUM_COLUMNS = 9;

using GffReader = io::CSVReader<GFF_NUM_COLUMNS, io::trim_chars<' '>, io::no_quote_escape<'\t'>, io::throw_on_overflow,
  io::single_and_empty_line_comment<'#', ' '> >;

using AttribMap = std::map<std::string, std::string>;

constexpr const auto GFF_COLUMNS_REQUIRED = {
  "feature",
  "start",
  "end",
  "strand",
  "frame",
  "attribute",
};

constexpr const auto GFF_ATTRIB_KEY_GENE_NAME = "gene_name";

class ErrorGffParserMissingColumn : public std::runtime_error {
public:
  explicit ErrorGffParserMissingColumn(const std::string& colName)
      : std::runtime_error(fmt::format("GFF parser: column \"{:s}\" is missing", colName)) {}
};

class ErrorGffParserInvalidAttribute : public std::runtime_error {
public:
  explicit ErrorGffParserInvalidAttribute(const std::string& attrib)
      : std::runtime_error(fmt::format(
          "GFF parser: unable to parse attribute: expected a key-value pair format is `key \"value\"`, got `{:s}`,",
          attrib)) {}
};


class ErrorGffParserMissingGeneName : public std::runtime_error {
public:
  explicit ErrorGffParserMissingGeneName()
      : std::runtime_error(fmt::format("GFF parser: expected attribute `${:s}` is missing", GFF_ATTRIB_KEY_GENE_NAME)) {
  }
};

class ErrorGffParserGenNameEmpty : public std::runtime_error {
public:
  explicit ErrorGffParserGenNameEmpty()
      : std::runtime_error(fmt::format("GFF parser: empty gene names are not allowed")) {}
};

class ErrorGffParserGeneRangeInvalid : public std::runtime_error {
public:
  explicit ErrorGffParserGeneRangeInvalid(const std::string& geneName, int start, int end)
      : std::runtime_error(
          fmt::format("GFF parser: in gene \"{:s}\": range is invalid: ({:d}, {:d})", geneName, start, end)) {}
};

class ErrorGffParserGeneFrameInvalid : public std::runtime_error {
public:
  explicit ErrorGffParserGeneFrameInvalid(const std::string& geneName, int frame)
      : std::runtime_error(fmt::format(
          R"(GFF parser: in gene "{:s}": frame "{:d}" is invalid, expected an integer 0, 1, or 2)", geneName, frame)) {}
};


class ErrorGffParserGeneStrandInvalid : public std::runtime_error {
public:
  explicit ErrorGffParserGeneStrandInvalid(const std::string& geneName, const std::string& strand)
      : std::runtime_error(fmt::format(
          R"(GFF parser: in gene "{:s}": gene strand is invalid: "{:s}", expected "+" or "-")", geneName, strand)) {}
};


std::pair<std::string, std::string> parseAttribute(const std ::string& keyValStr) {
  std::vector<std::string> keyVal;
  boost::split(keyVal, keyValStr, boost::is_any_of(" "));

  if (keyVal.size() != 2) {
    throw ErrorGffParserInvalidAttribute(keyValStr);
  }

  auto& key = keyVal[0];
  key = boost::remove_erase_if(key, boost::is_any_of(R"('" )"));

  auto& val = keyVal[1];
  val = boost::remove_erase_if(val, boost::is_any_of(R"('" )"));

  if (key.empty() || val.empty()) {
    throw ErrorGffParserInvalidAttribute(keyValStr);
  }

  // GFF allows for empty values if they are filled with dots. We don't need these dots.
  if (val == ".") {
    val = "";
  }

  return std::make_pair(key, val);
}

AttribMap parseAttributes(const std ::string& attribsStr) {
  std::vector<std::string> keyValStrs;
  boost::split(keyValStrs, attribsStr, boost::is_any_of(";"));

  AttribMap attribMap;
  std::transform(keyValStrs.begin(), keyValStrs.end(), std::inserter(attribMap, attribMap.end()), parseAttribute);

  return attribMap;
}


std::string getGeneName(const AttribMap& attribMap) {
  const auto& found = attribMap.find(GFF_ATTRIB_KEY_GENE_NAME);
  if (found == attribMap.end()) {
    throw ErrorGffParserMissingGeneName();
  }
  return found->second;
}

void validateGene(const Gene& gene) {
  if (gene.geneName.empty()) {
    throw ErrorGffParserGenNameEmpty();
  }

  if (gene.start < 1 || gene.start > gene.end) {
    throw ErrorGffParserGeneRangeInvalid(gene.geneName, gene.start, gene.end);
  }

  if (gene.frame < 0 || gene.frame > 2) {
    throw ErrorGffParserGeneFrameInvalid(gene.geneName, gene.frame);
  }

  if (gene.strand != "+" && gene.strand != "-") {
    throw ErrorGffParserGeneStrandInvalid(gene.geneName, gene.strand);
  }
}

GeneMap parseGeneMapGff(std::istream& is) {
  GffReader reader("file", is);

  reader.set_header(
    /* 1 */ "seqname",
    /* 2 */ "source",
    /* 3 */ "feature",
    /* 4 */ "start",
    /* 5 */ "end",
    /* 6 */ "score",
    /* 7 */ "strand",
    /* 8 */ "frame",
    /* 9 */ "attribute");

  std::for_each(std::begin(GFF_COLUMNS_REQUIRED), std::end(GFF_COLUMNS_REQUIRED), [&reader](const auto& colName) {
    if (!reader.has_column(colName)) {
      throw ErrorGffParserMissingColumn(colName);
    }
  });

  /* 1 */ [[maybe_unused]] std::string seqname;
  /* 2 */ [[maybe_unused]] std::string source;
  /* 3 */ std::string feature;
  /* 4 */ int start;// NOLINT(cppcoreguidelines-init-variables)
  /* 5 */ int end;  // NOLINT(cppcoreguidelines-init-variables)
  /* 6 */ [[maybe_unused]] std::string score;
  /* 7 */ std::string strand;
  /* 8 */ int frame;// NOLINT(cppcoreguidelines-init-variables)
  /* 9 */ [[maybe_unused]] std::string attribute;

  GeneMap geneMap;
  while (reader.read_row(seqname, source, feature, start, end, score, strand, frame, attribute)) {
    if (feature != "gene") {
      continue;
    }

    const auto attribMap = parseAttributes(attribute);
    const auto& geneName = getGeneName(attribMap);

    const Gene gene = {
      .geneName = geneName,
      .start = start,
      .end = end,
      .strand = strand,
      .frame = frame,
    };

    validateGene(gene);

    geneMap.emplace(geneName, gene);
  }

  return geneMap;
}
