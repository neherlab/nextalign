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

find_package(boost_base 1.68.0 REQUIRED)
find_package(fast-cpp-csv-parser 20191004 REQUIRED)
find_package(fmt 7.1.0 REQUIRED)
find_package(Microsoft.GSL 3.1.0 REQUIRED)

add_library(${PROJECT_NAME} STATIC
  include/nextalign/geneUtils.h
  include/nextalign/nextalign.h
  include/nextalign/parseGb.h
  include/nextalign/parseGeneMapGff.h
  include/nextalign/parseSequences.h
  include/nextalign/types.h
  src/alignPairwise.cpp
  src/alignPairwise.h
  src/geneUtils.cpp
  src/matchNuc.cpp
  src/matchNuc.h
  src/nextalign.cpp
  src/parseGb.cpp
  src/parseGeneMapGff.cpp
  src/parseSequences.cpp
  src/safeCast.h
  src/seqan2.cpp
  src/seqan2.h
  src/seqan2nowarn.h
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
  "${CMAKE_CURRENT_SOURCE_DIR}/../../3rdparty/seqan/include"
  )

target_link_libraries(${PROJECT_NAME}
  PUBLIC

  PRIVATE
  fast-cpp-csv-parser::fast-cpp-csv-parser
  boost_base::boost_base
  fmt::fmt
  Microsoft.GSL::GSL
  )

set(${PROJECT_NAME}_BUILD_TESTS ON CACHE BOOL "Build unit tests")
if (${${PROJECT_NAME}_BUILD_TESTS})
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif ()
