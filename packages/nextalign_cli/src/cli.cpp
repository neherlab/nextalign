#include <fmt/format.h>
#include <nextalign/nextalign.h>
#include <nextalign/parseGeneMapGff.h>
#include <nextalign/parseSequences.h>
#include <nextalign/types.h>

#include <boost/algorithm/string.hpp>
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
  std::string genes;
  std::string output;
};


template<typename Result>
auto getParamRequired(
  const cxxopts::Options &cxxOpts, const cxxopts::ParseResult &cxxOptsParsed, const std::string &name) -> Result {
  if (!cxxOptsParsed.count(name)) {
    fmt::print(stderr, "Error: argument `--{:s}` is required\n\n", name);
    fmt::print(stderr, "{:s}\n", cxxOpts.help());
    std::exit(1);
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
       "(required) List of genes to account for during alignment",
       cxxopts::value<std::string>(),
       "GENES"
    )

    (
      "o,output",
      "(required) Path to output aligned sequence in FASTA format",
      cxxopts::value<std::string>(),
      "OUTPUT"
    );
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
  const auto genes = getParamRequired<std::string>(cxxOpts, cxxOptsParsed, "genes");
  const auto output = getParamRequired<std::string>(cxxOpts, cxxOptsParsed, "output");

  return {
    jobs,
    sequences,
    reference,
    genemap,
    genes,
    output,
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

std::set<std::string> parseGenes(const CliParams &cliParams) {
  std::set<std::string> genes;
  boost::algorithm::split(genes, cliParams.genes, boost::is_any_of(","));
  return genes;
}

void validateGenes(const std::set<std::string> &genes, const GeneMap &geneMap) {
  if (genes.empty()) {
    fmt::print("Warning: no target genes provided for codon-aware alignment. Alignment quality is degraded.\n");
  }

  for (const auto &gene : genes) {
    const auto &it = geneMap.find(gene);
    if (it == geneMap.end()) {
      fmt::print(stderr, "Error: gene \"{}\" is not in gene map", gene);
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
  fmt::format_to(buf, "{:>12s}=\"{:<s}\"\n", "--genes", cliParams.genes);
  fmt::format_to(buf, "{:>12s}=\"{:<s}\"\n", "--output", cliParams.output);
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


int main(int argc, char *argv[]) {
  try {
    const auto cliParams = parseCommandLine(argc, argv);
    fmt::print(stdout, formatCliParams(cliParams));

    NextalignOptions options;

    const auto [refName, ref] = parseRefFastaFile(cliParams.reference);
    fmt::print(stdout, formatRef(refName, ref));

    const auto geneMap = parseGeneMapGffFile(cliParams.genemap);
    options.genes = parseGenes(cliParams);
    validateGenes(options.genes, geneMap);

    fmt::print(stdout, formatGeneMap(geneMap, options.genes));

    std::ifstream fastaFile(cliParams.sequences);
    const auto fastaStream = makeFastaStream(fastaFile);
    if (!fastaFile.good()) {
      fmt::print(stderr, "Error: unable to read \"{:s}\"\n", cliParams.sequences);
      std::exit(1);
    }


    constexpr const auto TABLE_WIDTH = 71;
    fmt::print(stdout, "\nSequences:\n");
    fmt::print(stdout, "{:s}\n", std::string(TABLE_WIDTH, '-'));
    fmt::print(stdout, "| {:5s} | {:40s} | {:16s} |\n", "Index", "Seq. name", "Align. score");
    fmt::print(stdout, "{:s}\n", std::string(TABLE_WIDTH, '-'));
    int i = 0;
    while (fastaStream->good()) {
      const auto entry = fastaStream->next();
      fmt::print(stdout, "| {:5d} | {:<40s} | ", i, entry.first);

      try {
        const auto &alignment = nextalign(entry.second, ref, geneMap, options);
        fmt::print(stdout, "{:>16d} |\n", alignment.alignmentScore);
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
