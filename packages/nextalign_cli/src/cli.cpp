#include <fmt/format.h>
#include <nextalign/nextalign.h>
#include <tbb/global_control.h>
#include <tbb/pipeline.h>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string/split.hpp>
#include <cxxopts.hpp>
#include <filesystem>
#include <fstream>


struct CliParams {
  int jobs;
  std::string sequences;
  std::string reference;
  std::string genemap;
  std::optional<std::string> genes;
  std::optional<std::string> outputDir;
  std::optional<std::string> outputBasename;
  std::optional<std::string> outputFasta;
  std::optional<std::string> outputInsertions;
};

struct Paths {
  std::filesystem::path outputFasta;
  std::filesystem::path outputInsertions;
};

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
      "(optional) Number of CPU threads used by the algorithm. If not specified or non-positive, will use all available threads",
      cxxopts::value<int>()->default_value(std::to_string(0)),
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
      "d,output-dir",
      "(optional) Write output files to this directory. The base filename can be set using --output-basename flag. The paths can be overridden on a per-file basis using --output-* flags. If the required directory tree does not exist, it will be created.",
      cxxopts::value<std::string>(),
      "OUTPUT"
    )

    (
      "n,output-basename",
      "(optional) Sets the base filename to use for output files. To be used together with --output-dir flag. By default uses the filename of the sequences file (provided with --sequences). The paths can be overridden on a per-file basis using --output-* flags.",
      cxxopts::value<std::string>(),
      "OUTPUT"
    )

    (
      "o,output-fasta",
      "(required) Path to output aligned sequences in FASTA format (overrides paths given with --output-dir and --output-basename). If the required directory tree does not exist, it will be created.",
      cxxopts::value<std::string>(),
      "OUTPUT"
    )

    (
      "I,output-insertions",
      "(optional) Path to output stripped insertions data in CSV format (overrides paths given with --output-dir and --output-basename). If the required directory tree does not exist, it will be created.",
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
  const auto outputDir = getParamOptional<std::string>(cxxOpts, cxxOptsParsed, "output-dir");
  const auto outputBasename = getParamOptional<std::string>(cxxOpts, cxxOptsParsed, "output-basename");
  const auto outputFasta = getParamOptional<std::string>(cxxOpts, cxxOptsParsed, "output-fasta");
  const auto outputInsertions = getParamOptional<std::string>(cxxOpts, cxxOptsParsed, "output-insertions");

  return {
    jobs,
    sequences,
    reference,
    genemap,
    genes,
    outputDir,
    outputBasename,
    outputFasta,
    outputInsertions,
  };
}

AlgorithmInput parseRefFastaFile(const std::string &filename) {
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
  fmt::format_to(buf, "{:>20s}=\"{:<d}\"\n", "--jobs", cliParams.jobs);
  fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "--sequences", cliParams.sequences);
  fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "--reference", cliParams.reference);
  fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "--genemap", cliParams.genemap);

  if (cliParams.genes) {
    fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "--genes", *cliParams.genes);
  }

  if (cliParams.outputDir) {
    fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "--output-dir", *cliParams.outputDir);
  }

  if (cliParams.outputBasename) {
    fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "--output-basename", *cliParams.outputBasename);
  }

  if (cliParams.outputFasta) {
    fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "--output-fasta", *cliParams.outputFasta);
  }

  if (cliParams.outputInsertions) {
    fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "--output-insertions", *cliParams.outputInsertions);
  }

  return fmt::to_string(buf);
}

Paths getPaths(const CliParams &cliParams) {
  std::filesystem::path sequencesPath = cliParams.sequences;

  auto outDir = std::filesystem::canonical(std::filesystem::current_path());
  if (cliParams.outputDir) {
    outDir = *cliParams.outputDir;
  }

  if (!outDir.is_absolute()) {
    outDir = std::filesystem::current_path() / outDir;
  }

  fmt::print("OUT PATH: \"{:<s}\"\n", outDir.string());

  auto baseName = sequencesPath.stem();
  if (cliParams.outputBasename) {
    baseName = *cliParams.outputBasename;
  }

  auto outputFasta = outDir / baseName;
  outputFasta += ".aligned.fasta";
  if (cliParams.outputFasta) {
    outputFasta = *cliParams.outputFasta;
  }

  auto outputInsertions = outDir / baseName;
  outputInsertions += ".insertions.csv";
  if (cliParams.outputInsertions) {
    outputInsertions = *cliParams.outputInsertions;
  }

  return {.outputFasta = outputFasta, .outputInsertions = outputInsertions};
}

std::string formatPaths(const Paths &paths) {
  fmt::memory_buffer buf;
  fmt::format_to(buf, "\nOutput files:\n");
  fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "Aligned sequences", paths.outputFasta.string());
  fmt::format_to(buf, "{:>20s}=\"{:<s}\"\n", "Stripped insertions", paths.outputInsertions.string());
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

template<typename Letter>
inline constexpr auto is_letter(Letter left) {
  return [left](Letter right) { return left == right; };
}

/**
 * Runs nextalign algorithm in a parallel pipeline
 */
void run(
  /* in  */ int parallelism,
  /* in  */ const CliParams &cliParams,
  /* out */ std::unique_ptr<FastaStream> &inputFastaStream,
  /* in  */ const std::string &refStr,
  /* in  */ const GeneMap &geneMap,
  /* in  */ const NextalignOptions &options,
  /* out */ std::ostream &outputFastaStream,
  /* out */ std::ostream &outputInsertionsStream) {
  tbb::task_group_context context;

  const auto ref = toNucleotideSequence(refStr);

  /** Input filter is a serial input filter function, which accepts an input stream,
   * reads and parses the contents of it, and returns parsed sequences */
  const auto inputFilter = tbb::make_filter<void, AlgorithmInput>(tbb::filter::serial_in_order,//
    [&inputFastaStream](tbb::flow_control &fc) -> AlgorithmInput {
      if (!inputFastaStream->good()) {
        fc.stop();
        return {};
      }

      return inputFastaStream->next();
    });


  /** A set of parallel transform filter functions, each accepts a parsed sequence from the input filter,
   * runs nextalign algorithm sequentially and returning its result.
   * The number of filters is determined by the `--jobs` CLI argument */
  const auto transformFilters = tbb::make_filter<AlgorithmInput, AlgorithmOutput>(tbb::filter::parallel,//
    [&ref, &geneMap, &options](const AlgorithmInput &input) -> AlgorithmOutput {
      try {
        auto query = toNucleotideSequence(input.seq);
        boost::trim_if(query, is_letter(Nucleotide::N));

        const auto result = nextalign(query, ref, geneMap, options);
        return {.index = input.index, .seqName = input.seqName, .hasError = false, .result = result, .error = nullptr};
      } catch (const std::exception &e) {
        const auto &error = std::current_exception();
        return {.index = input.index, .seqName = input.seqName, .hasError = true, .result = {}, .error = error};
      }
    });


  /** Output filter is a serial ordered filter function which accepts the results from transform filters,
   * one at a time, displays and writes them to output streams */
  const auto outputFilter = tbb::make_filter<AlgorithmOutput, void>(tbb::filter::serial_in_order,//
    [&cliParams, &outputFastaStream, &outputInsertionsStream](const AlgorithmOutput &output) {
      const auto index = output.index;
      const auto &seqName = output.seqName;

      const auto &error = output.error;
      if (error) {
        try {
          std::rethrow_exception(error);
        } catch (const std::exception &e) {
          fmt::print(stdout, "{:s}:{:s}\n", seqName, e.what());
          return;
        }
      }

      const auto &query = output.result.query;
      const auto &alignmentScore = output.result.alignmentScore;
      const auto &insertions = output.result.insertions;
      fmt::print(stdout, "| {:5d} | {:<40s} | {:>16d} | {:12d} | \n",//
        index, seqName, alignmentScore, insertions.size());

      outputFastaStream << fmt::format(">{:s}\n{:s}\n", seqName, query);

      outputInsertionsStream << fmt::format("\"{:s}\",\"{:s}\"\n", seqName, formatInsertions(insertions));
    });

  try {
    tbb::parallel_pipeline(parallelism, inputFilter & transformFilters & outputFilter, context);
  } catch (const std::exception &e) {
    fmt::print(stdout, "{:s}\n", e.what());
  }
}


int main(int argc, char *argv[]) {
  try {
    const auto cliParams = parseCommandLine(argc, argv);
    fmt::print(stdout, formatCliParams(cliParams));

    const auto paths = getPaths(cliParams);
    fmt::print(stdout, formatPaths(paths));

    std::filesystem::create_directories(paths.outputFasta.parent_path());
    std::filesystem::create_directories(paths.outputInsertions.parent_path());

    NextalignOptions options;

    const auto refInput = parseRefFastaFile(cliParams.reference);
    const auto &refName = refInput.seqName;
    const auto &ref = refInput.seq;
    fmt::print(stdout, formatRef(refName, ref));

    const auto geneMap = parseGeneMapGffFile(cliParams.genemap);
    options.genes = parseGenes(cliParams, geneMap);
    validateGenes(options.genes, geneMap);

    fmt::print(stdout, formatGeneMap(geneMap, options.genes));

    std::ifstream fastaFile(cliParams.sequences);
    auto fastaStream = makeFastaStream(fastaFile);
    if (!fastaFile.good()) {
      fmt::print(stderr, "Error: unable to read \"{:s}\"\n", cliParams.sequences);
      std::exit(1);
    }

    std::ofstream outputFastaFile(paths.outputFasta);
    if (!outputFastaFile.good()) {
      fmt::print(stderr, "Error: unable to write \"{:s}\"\n", paths.outputFasta.string());
      std::exit(1);
    }


    std::ofstream outputInsertionsFile(paths.outputInsertions);
    if (!outputInsertionsFile.good()) {
      fmt::print(stderr, "Error: unable to write \"{:s}\"\n", paths.outputInsertions.string());
      std::exit(1);
    }
    outputInsertionsFile << "seqName,insertions\n";

    int parallelism = -1;
    if (cliParams.jobs > 0) {
      tbb::global_control globalControl{
        tbb::global_control::max_allowed_parallelism, static_cast<size_t>(cliParams.jobs)};
      parallelism = cliParams.jobs;
    }

    fmt::print("\nParallelism: {:d}\n", parallelism);

    constexpr const auto TABLE_WIDTH = 86;
    fmt::print(stdout, "\nSequences:\n");
    fmt::print(stdout, "{:s}\n", std::string(TABLE_WIDTH, '-'));
    fmt::print(stdout, "| {:5s} | {:40s} | {:16s} | {:12s} |\n", "Index", "Seq. name", "Align. score", "Insertions");
    fmt::print(stdout, "{:s}\n", std::string(TABLE_WIDTH, '-'));


    try {
      run(parallelism, cliParams, fastaStream, ref, geneMap, options, outputFastaFile, outputInsertionsFile);
    } catch (const std::exception &e) {
      fmt::print(stdout, "{:>16s} |\n", e.what());
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
