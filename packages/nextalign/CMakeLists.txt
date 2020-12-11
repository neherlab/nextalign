cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0069 NEW)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/packages/cmake" ${CMAKE_MODULE_PATH})
include(CFlags)
include(Quiet)
include(Sanitizers)

project(nextalign)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_MESSAGE LAZY)

find_package(Microsoft.GSL 3.1.0 REQUIRED)
find_package(boost_base 1.68.0 REQUIRED)
find_package(fast-cpp-csv-parser 20191004 REQUIRED)
find_package(fmt 7.1.0 REQUIRED)

add_library(${PROJECT_NAME} STATIC
  include/nextalign/geneUtils.h
  include/nextalign/nextalign.h
  include/nextalign/parseGb.h
  include/nextalign/parseGeneMapGff.h
  include/nextalign/parseSequences.h
  include/nextalign/types.h
  src/alignCodon.cpp
  src/alignCodon.h
  src/alignPairwise.cpp
  src/alignPairwise.h
  src/decode.cpp
  src/decode.h
  src/geneUtils.cpp
  src/matchNuc.cpp
  src/matchNuc.h
  src/nextalign.cpp
  src/parseGb.cpp
  src/parseGeneMapGff.cpp
  src/parseSequences.cpp
  src/safeCast.h
  src/translate.cpp
  src/translate.h
  src/utils/config.h
  src/utils/contract.h
  src/utils/debugbreak.h
  src/utils/ieee754_comparison.h
  src/utils/macro_overload.h
  )

set_property(TARGET ${PROJECT_NAME} PROPERTY C_STANDARD 11)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

target_compile_definitions(${PROJECT_NAME}
  PRIVATE -DPROJECT_NAME="${PROJECT_NAME}"
  )

target_include_directories(${PROJECT_NAME}
  PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src"

  PUBLIC
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  )

target_include_directories(${PROJECT_NAME} SYSTEM
  PRIVATE
  "${CMAKE_SOURCE_DIR}/3rdparty/frozen/include"
  )

target_link_libraries(${PROJECT_NAME}
  PUBLIC

  PRIVATE
  Microsoft.GSL::GSL
  boost_base::boost_base
  fast-cpp-csv-parser::fast-cpp-csv-parser
  fmt::fmt
  )

set(NEXTALIGN_BUILD_TESTS ON CACHE BOOL "Build unit tests")
if (${NEXTALIGN_BUILD_TESTS})
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif ()

set(NEXTALIGN_BUILD_BENCHMARKS OFF CACHE BOOL "Build benchmarks")
if (${NEXTALIGN_BUILD_BENCHMARKS})
  add_subdirectory(benchmarks)
endif ()
