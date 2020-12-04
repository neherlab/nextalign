<h1 id="nextclade" align="center">
Nextalign
</h1>

<h4 id="nextclade" align="center">
Sequence alignment
</h1>

## Developer's guide

### Quick start

1. Install and configure required dependencies

    - cmake >= 3.10

    - conan package manager (or Ubuntu Linux there is an installation script included in `./tools/install-conan`)

    - (optional, but recommended) nodemon (requires Node.js)

   ```bash
   npm install --global nodemon
   ```


2. Clone and run

   ```bash
   git clone --recursive https://github.com/neherlab/nextalign
   cd nextalign
   make dev

   ```
   (note the `--recursive` flag for git - this repositoy contains git submodules)

   This will:

    - install or update conan packages
    - run cmake and generate makefiles
    - build the project and tests
    - run static analysis on source files
    - run tests
    - watch source files and rebuild on changes

   If you don't want to install Node.js and nodemon, or don't want the automatic rebuild, you can use `make dev-nowatch`
   instead of `make dev`.

### Benchmarks

A set of benchmarks is located in [`packages/nextalign/benchmarks`](https://github.com/neherlab/nextalign/tree/master/packages/nextalign/benchmarks). We are
using [Google Benchmark](https://github.com/google/benchmark) framework. Read the important [Runtime and Reporting Considerations](https://github.com/google/benchmark#runtime-and-reporting-considerations).

For the most accurate results, [disable CPU frequence scaling](https://github.com/google/benchmark#disabling-cpu-frequency-scaling) for the time of
your benchmarking session. (More info: [[kernel](https://www.kernel.org/doc/html/v4.15/admin-guide/pm/cpufreq.html)]
, [[arch](https://wiki.archlinux.org/index.php/CPU_frequency_scaling)]
, [[debian](https://wiki.debian.org/CpuFrequencyScaling)])

Run benchmarks with

```bash
make benchmarks
```

This will install dependencies, build the library and benchmarks in "Release" mode and will run the benchmarks.
Benchmarks will rerun on code changes.

Or run the `scripts/benchmarks.sh` directly (no hot reloading).

For better debugging, you can also build & run in "Debug" mode and under GDB with:

```bash
CMAKE_BUILD_TYPE=Debug USE_GDB=1 make benchmarks
```

You can pass parameters to the benchmark executable with either of:

```bash
BENCHMARK_OPTIONS='--help' make benchmarks
scripts/benchmarks.sh --help
```

For example, you can filter the benchmarks by name: to run only the benchmarks containing the word "Average":

```bash
BENCHMARK_OPTIONS='--benchmark_filter=.*Average' make benchmarks
```

The results are also saved to the file `nextalign_benchmarks.json`.
You can compare multiple results using the [compare.py](https://github.com/google/benchmark/tree/master/tools) tool in the Google Benchmark's repository. For more information refer to [Benchmark Tools](https://github.com/google/benchmark/blob/master/docs/tools.md) documentation.


### Tests

Test are run as a part of the main development script (`make dev`). We are using [Google Test](https://github.com/google/googletest/blob/master/googletest/docs/primer.md)
See [Google Test documentation](https://github.com/google/googletest/blob/master/googletest/docs/primer.md) and [Google Mock documentation](https://github.com/google/googletest/blob/master/googlemock/README.md) for more details.

## License

<a target="_blank" rel="noopener noreferrer" href="LICENSE" alt="License file">MIT License</a>
