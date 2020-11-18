#include <fmt/format.h>
#include <nextalign/nextalign.h>
#include <nextalign/parseGb.h>
#include <nextalign/types.h>

#include <cxxopts.hpp>

// TODO(ivan-aksamentov): detect number of cores
const int numCores = 4;


using CliOptions = cxxopts::ParseResult;

auto parseCommandLine(int argc, char *argv[]) {// NOLINT(cppcoreguidelines-avoid-c-arrays)
  cxxopts::Options options("nextalign", "Nextalign: sequence alignment");

  // clang-format off
  options.add_options()
    (
      "h,help",
      "Show this help"
    )

    (
      "j,jobs",
      "Number of CPU threads used by the algorithm. If not specified, using number of available logical CPU cores",
      cxxopts::value<int>()->default_value(std::to_string(numCores)),
      "JOBS"
    )

    (
      "i,sequences",
      "Path to a FASTA or file with input sequences",
      cxxopts::value<std::string>(),
      "SEQS"
    )

    (
      "r,reference",
       R"(Path to a GB file containing reference sequence and gene map)",
       cxxopts::value<std::string>(),
       "REF"
    )

    (
      "m,genemap",
       R"(Path to a JSON file containing custom gene map)",
       cxxopts::value<std::string>(),
       "GENEMAP"
    )

    (
      "g,genes",
       R"(List of genes to account for during alignment)",
       cxxopts::value<std::string>(),
       "GENES"
    )

    (
      "o,output",
      "(optional) Path to output aligned sequence ins FASTA format",
      cxxopts::value<std::string>(),
      "OUTPUT"
    );
  // clang-format on

  auto result = options.parse(argc, argv);

  if (result.count("help") > 0) {
    std::cout << options.help() << std::endl;
    std::exit(0);
  }

  if (result.count("input-fasta") == 0) {
    std::cerr << "Error: input-fasta argument is required" << std::endl;
    std::cerr << options.help() << std::endl;
    std::exit(1);
  }

  return result;
}

class DeveloperError : public std::runtime_error {
public:
  explicit DeveloperError(const std::string &message) : std::runtime_error(message) {}
};

template<typename Result>
auto getOptionRequired(const CliOptions &params, const std::string &name) -> Result {
  if (!params.count(name)) {
    throw DeveloperError(
      fmt::format("Developer Error: `--{:s}` argument is required but is not present. This may mean that "
                  "the missing required argument has escaped the argument validation. This needs to be fixed.",
        name));
  }

  return params[name].as<Result>();
}


struct CliParams {
  int numThreads;
  std::string inputFasta;
  std::string inputRootSeq;
  std::string inputGeneMap;
  std::string inputGenes;
};

CliParams validateCliParams(const CliOptions &options) {
  const auto numThreads = getOptionRequired<int>(options, "jobs");
  const auto inputFasta = getOptionRequired<std::string>(options, "sequences");
  const auto inputRootSeq = getOptionRequired<std::string>(options, "root-seq");
  const auto inputGeneMap = getOptionRequired<std::string>(options, "gene-map");
  const auto inputGenes = getOptionRequired<std::string>(options, "genes");
  const auto outputFasta = getOptionRequired<std::string>(options, "output");
  return {numThreads, inputFasta, inputRootSeq, inputGeneMap, outputFasta};
}


int main(int argc, char *argv[]) {
  try {
    const auto args = parseCommandLine(argc, argv);
    const auto params = validateCliParams(args);

    std::cout << "numThreads   : " << params.numThreads << std::endl;
    std::cout << "inputFasta   : " << params.inputFasta << std::endl;
    std::cout << "inputRootSeq : " << params.inputRootSeq << std::endl;
    std::cout << "inputGeneMap : " << params.inputGeneMap << std::endl;
    std::cout << "inputGenes   : " << params.inputGenes << std::endl;


    const std::string gbContent = "TODO";
    const NextalignOptions options = {};

    // Parse and prepare reference sequence and genemap
    const auto [ref, geneMap] = parseGb(gbContent);

    // for each streamed sequence
    for (;;) {
      std::string query("TODO");
      nextalign(query, ref, geneMap, options);
    }

  } catch (const cxxopts::OptionSpecException &e) {
    std::cerr << e.what() << std::endl;
  } catch (const cxxopts::OptionParseException &e) {
    std::cerr << e.what() << std::endl;
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
