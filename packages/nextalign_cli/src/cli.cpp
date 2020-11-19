#include <fmt/format.h>
#include <nextalign/nextalign.h>
#include <nextalign/parseGb.h>
#include <nextalign/parseSequences.h>
#include <nextalign/types.h>

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


int main(int argc, char *argv[]) {
  try {
    const auto cliParams = parseCommandLine(argc, argv);

    std::cout << "Parameters" << std::endl;
    std::cout << "  jobs     : " << cliParams.jobs << std::endl;
    std::cout << "  sequences: " << cliParams.sequences << std::endl;
    std::cout << "  reference: " << cliParams.reference << std::endl;
    std::cout << "  genemap  : " << cliParams.genemap << std::endl;
    std::cout << "  genes    : " << cliParams.genes << std::endl;
    std::cout << "  output   : " << cliParams.output << std::endl;
    std::cout << std::endl;

    // const std::string gbContent = "TODO";
    const NextalignOptions options = {};

    // Parse and prepare reference sequence and genemap
    // const auto [ref, geneMap] = parseGb(gbContent);

    std::ifstream file(cliParams.sequences);
    const auto fastaStream = makeFastaStream(file);
    if (!file.good()) {
      fmt::print(stderr, "Error: unable to read {:s}\n", cliParams.sequences);
      std::exit(1);
    }

    while (fastaStream->good()) {
      const auto entry = fastaStream->next();
      fmt::print(stdout, "{:s}\n", entry.first);
    }

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
