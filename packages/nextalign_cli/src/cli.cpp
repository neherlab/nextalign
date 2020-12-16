#include <fmt/format.h>
#include <nextalign/nextalign.h>
#include <nextalign/parseGeneMapGff.h>
#include <nextalign/parseSequences.h>
#include <nextalign/types.h>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cxxopts.hpp>
#include <fstream>

// TODO(ivan-aksamentov): detect number of cores
const int numCores = 4;


struct CliParams {
  int jobs;
  std::string sequences;
  std::string reference;
  std::string genemap;
  std::optional<std::string> genes;
  std::string output;
  std::optional<std::string> outputInsertions;
};


std::string formatInsertions(const std::vector<Insertion> &insertions);
template<typename Result>
Result getParamRequired(
  const cxxopts::Options &cxxOpts, const cxxopts::ParseResult &cxxOptsParsed, const std::string &name) {
  if (!cxxOptsParsed.count(name)) {
    fmt::print(stderr, "Error: argument `--{:s}` is required\n\n", name);
    fmt::print(stderr, "{:s}\n", cxxOpts.help());
    std::exit(1);
  }

  return cxxOptsParsed[name].as<Result>();
}

template<typename Result>
std::optional<Result> getParamOptional(
  const cxxopts::Options &cxxOpts, const cxxopts::ParseResult &cxxOptsParsed, const std::string &name) {
  if (!cxxOptsParsed.count(name)) {
    return {};
  }

  return cxxOptsParsed[name].as<Result>();
}

template<typename Result>
auto getParamRequiredDefaulted([[maybe_unused]] const cxxopts::Options &cxxOpts,
  const cxxopts::ParseResult &cxxOptsParsed, const std::string &name) -> Result {
  return cxxOptsParsed[name].as<Result>();
}


CliParams parseCommandLine(int argc, char *argv[]) {// NOLINT(cppcoreguidelines-avoid-c-arrays)
  cxxopts::Options cxxOpts("nextalign", "Nextalign: sequence alignment\n");

  // clang-format off
  cxxOpts.add_options()
    (
      "h,help",
      "Show this help"
    )

    (
      "j,jobs",
      "(optional) Number of CPU threads used by the algorithm. If not specified, using number of available logical CPU cores",
      cxxopts::value<int>()->default_value(std::to_string(numCores)),
      "JOBS"
    )

    (
      "i,sequences",
      "(required) Path to a FASTA or file with input sequences",
      cxxopts::value<std::string>(),
      "SEQS"
    )

    (
      "r,reference",
       "(required) Path to a GB file containing reference sequence and gene map",
       cxxopts::value<std::string>(),
       "REF"
    )

    (
      "m,genemap",
       "(required) Path to a JSON file containing custom gene map",
       cxxopts::value<std::string>(),
       "GENEMAP"
    )

    (
      "g,genes",
       "(optional) List of genes to account for during alignment. If not supplied or empty, all the genes from the gene map are used",
       cxxopts::value<std::string>(),
       "GENES"
    )

    (
      "o,output",
      "(required) Path to output aligned sequences in FASTA format",
      cxxopts::value<std::string>(),
      "OUTPUT"
    )

    (
      "I,output-insertions",
      "(optional) Path to output stripped insertions data",
      cxxopts::value<std::string>(),
      "OUTPUT_INSERTIONS"
    )

;
  // clang-format on

  const auto cxxOptsParsed = cxxOpts.parse(argc, argv);

  if (cxxOptsParsed.count("help") > 0) {
    fmt::print(stdout, "{:s}\n", cxxOpts.help());
    std::exit(0);
  }

  const auto jobs = getParamRequiredDefaulted<int>(cxxOpts, cxxOptsParsed, "jobs");
  const auto sequences = getParamRequired<std::string>(cxxOpts, cxxOptsParsed, "sequences");
  const auto reference = getParamRequired<std::string>(cxxOpts, cxxOptsParsed, "reference");
  const auto genemap = getParamRequired<std::string>(cxxOpts, cxxOptsParsed, "genemap");
  const auto genes = getParamOptional<std::string>(cxxOpts, cxxOptsParsed, "genes");
  const auto output = getParamRequired<std::string>(cxxOpts, cxxOptsParsed, "output");
  const auto outputInsertions = getParamOptional<std::string>(cxxOpts, cxxOptsParsed, "output-insertions");

  return {
    jobs,
    sequences,
    reference,
    genemap,
    genes,
    output,
    outputInsertions,
  };
}

std::pair<std::string, std::string> parseRefFastaFile(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.good()) {
    fmt::print(stderr, "Error: unable to read \"{:s}\"\n", filename);
    std::exit(1);
  }

  const auto refSeqs = parseSequences(file);
  if (refSeqs.size() != 1) {
    fmt::print(stderr, "Error: {:d} sequences found in reference sequence file, expected 1", refSeqs.size());
    std::exit(1);
  }

  return *(refSeqs.begin());
}

GeneMap parseGeneMapGffFile(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.good()) {
    fmt::print(stderr, "Error: unable to read \"{:s}\"\n", filename);
    std::exit(1);
  }

  auto geneMap = parseGeneMapGff(file);
  if (geneMap.empty()) {
    fmt::print(stderr, "Error: gene map is empty");
    std::exit(1);
  }

  return geneMap;
}

std::set<std::string> parseGenes(const CliParams &cliParams, const GeneMap &geneMap) {
  std::set<std::string> genes;

  if (cliParams.genes && !(cliParams.genes->empty())) {
    boost::algorithm::split(genes, *(cliParams.genes), boost::is_any_of(","));
    return genes;
  }

  // If `genes` is empty, return all gene names from the gene map
  std::transform(                                  //
    geneMap.begin(), geneMap.end(),                //
    std::inserter(genes, genes.begin()),           //
    [](const auto &entry) { return entry.first; });//

  return genes;
}

void validateGenes(const std::set<std::string> &genes, const GeneMap &geneMap) {
  for (const auto &gene : genes) {
    const auto &it = geneMap.find(gene);
    if (it == geneMap.end()) {
      fmt::print(stderr, "Error: gene \"{}\" is not in gene map\n", gene);
      std::exit(1);
    }
  }
}

std::string formatCliParams(const CliParams &cliParams) {
  fmt::memory_buffer buf;
  fmt::format_to(buf, "\nCLI Parameters:\n");
  fmt::format_to(buf, "{:>12s}=\"{:<d}\"\n", "--jobs", cliParams.jobs);
  fmt::format_to(buf, "{:>12s}=\"{:<s}\"\n", "--sequences", cliParams.sequences);
  fmt::format_to(buf, "{:>12s}=\"{:<s}\"\n", "--reference", cliParams.reference);
  fmt::format_to(buf, "{:>12s}=\"{:<s}\"\n", "--genemap", cliParams.genemap);

  if (cliParams.genes) {
    fmt::format_to(buf, "{:>12s}=\"{:<s}\"\n", "--genes", *cliParams.genes);
  }

  fmt::format_to(buf, "{:>12s}=\"{:<s}\"\n", "--output", cliParams.output);

  if (cliParams.outputInsertions) {
    fmt::format_to(buf, "{:>12s}=\"{:<s}\"\n", "--output-insertions", *cliParams.outputInsertions);
  }

  return fmt::to_string(buf);
}


std::string formatRef(const std::string &refName, const std::string &ref) {
  return fmt::format("\nReference:\n  name: \"{:s}\"\n  length: {:d}\n", refName, ref.size());
}

std::string formatGeneMap(const GeneMap &geneMap, const std::set<std::string> &genes) {
  constexpr const auto TABLE_WIDTH = 86;

  fmt::memory_buffer buf;
  fmt::format_to(buf, "\nGene map:\n");
  fmt::format_to(buf, "{:s}\n", std::string(TABLE_WIDTH, '-'));
  fmt::format_to(buf, "| {:8s} | {:16s} | {:8s} | {:8s} | {:8s} | {:8s} | {:8s} |\n", "Selected", "   Gene Name",
    "  Start", "  End", " Length", "  Frame", " Strand");
  fmt::format_to(buf, "{:s}\n", std::string(TABLE_WIDTH, '-'));
  for (const auto &[geneName, gene] : geneMap) {
    const auto selected = std::find(genes.cbegin(), genes.cend(), geneName) != genes.cend();
    const auto selectedStr = selected ? "  yes" : " ";
    fmt::format_to(buf, "| {:8s} | {:16s} | {:8d} | {:8d} | {:8d} | {:8d} | {:8s} |\n", selectedStr, geneName,
      gene.start + 1, gene.end + 1, gene.length, gene.frame + 1, gene.strand);
  }
  fmt::format_to(buf, "{:s}\n", std::string(TABLE_WIDTH, '-'));
  return fmt::to_string(buf);
}

std::string formatInsertions(const std::vector<Insertion> &insertions) {
  std::vector<std::string> insertionStrings;
  for (const auto &insertion : insertions) {
    insertionStrings.emplace_back(fmt::format("{:d}:{:s}", insertion.begin, insertion.seq));
  }

  return boost::algorithm::join(insertionStrings, ";");
}

int main(int argc, char *argv[]) {
  try {
    const auto cliParams = parseCommandLine(argc, argv);
    fmt::print(stdout, formatCliParams(cliParams));

    NextalignOptions options;

    const auto [refName, ref] = parseRefFastaFile(cliParams.reference);
    fmt::print(stdout, formatRef(refName, ref));

    const auto geneMap = parseGeneMapGffFile(cliParams.genemap);
    options.genes = parseGenes(cliParams, geneMap);
    validateGenes(options.genes, geneMap);

    fmt::print(stdout, formatGeneMap(geneMap, options.genes));

    std::ifstream fastaFile(cliParams.sequences);
    const auto fastaStream = makeFastaStream(fastaFile);
    if (!fastaFile.good()) {
      fmt::print(stderr, "Error: unable to read \"{:s}\"\n", cliParams.sequences);
      std::exit(1);
    }

    std::ofstream outputFastaFile(cliParams.output);
    if (!outputFastaFile.good()) {
      fmt::print(stderr, "Error: unable to write \"{:s}\"\n", cliParams.output);
      std::exit(1);
    }

    std::ofstream outputInsertionsFile;
    if (cliParams.outputInsertions) {
      outputInsertionsFile.open(*cliParams.outputInsertions);
      if (!outputInsertionsFile.good()) {
        fmt::print(stderr, "Error: unable to write \"{:s}\"\n", *cliParams.outputInsertions);
        std::exit(1);
      }

      outputInsertionsFile << "seqName,insertions\n";
    }

    constexpr const auto TABLE_WIDTH = 86;
    fmt::print(stdout, "\nSequences:\n");
    fmt::print(stdout, "{:s}\n", std::string(TABLE_WIDTH, '-'));
    fmt::print(stdout, "| {:5s} | {:40s} | {:16s} | {:12s} |\n", "Index", "Seq. name", "Align. score", "Insertions");
    fmt::print(stdout, "{:s}\n", std::string(TABLE_WIDTH, '-'));
    int i = 0;
    while (fastaStream->good()) {
      const auto &[seqName, seq] = fastaStream->next();
      fmt::print(stdout, "| {:5d} | {:<40s} | ", i, seqName);

      try {
        const auto &alignment = nextalign(seq, ref, geneMap, options);
        const auto &insertions = alignment.insertions;
        fmt::print(stdout, "{:>16d} | {:12d} | \n", alignment.alignmentScore, alignment.insertions.size());
        outputFastaFile << fmt::format(">{:s}\n{:s}\n", seqName, alignment.query);

        if (cliParams.outputInsertions) {
          outputInsertionsFile << fmt::format("\"{:s}\",\"{:s}\"\n", seqName, formatInsertions(insertions));
        }
      } catch (const std::exception &e) {
        fmt::print(stdout, "{:>16s} |\n", e.what());
      }

      ++i;
    }

    fmt::print(stdout, "{:s}\n", std::string(TABLE_WIDTH, '-'));

  } catch (const cxxopts::OptionSpecException &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::exit(1);
  } catch (const cxxopts::OptionParseException &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::exit(1);
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    std::exit(1);
  }
}
