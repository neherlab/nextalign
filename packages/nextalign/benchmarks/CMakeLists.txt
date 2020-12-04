cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0069 NEW)

project(nextalign_benchmarks)

find_package(benchmark REQUIRED)
find_package(fmt 7.1.0 REQUIRED)
find_package(range-v3 0.11.0 REQUIRED)

add_executable(${PROJECT_NAME}
  alignPairwise.benchmark.cpp
  )

set_property(TARGET ${PROJECT_NAME} PROPERTY C_STANDARD 11)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

target_link_libraries(${PROJECT_NAME}
  PRIVATE
  nextalign
  benchmark::benchmark
  fmt::fmt
  range-v3::range-v3
  )
