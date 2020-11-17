cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0069 NEW)

project(nextalign)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_MESSAGE LAZY)

find_package(boost_base 1.68.0 REQUIRED)
find_package(fmt 7.1.0 REQUIRED)

add_library(${PROJECT_NAME} STATIC
  src/alignPairwise.cpp
  src/alignPairwise.h
  src/parseSequences.cpp
  src/parseSequences.h
  )

set_property(TARGET ${PROJECT_NAME} PROPERTY C_STANDARD 11)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 14)

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
  boost_base::boost_base

  PRIVATE
  fmt::fmt
  )

set(${PROJECT_NAME}_BUILD_TESTS ON CACHE BOOL "Build unit tests")
if (${${PROJECT_NAME}_BUILD_TESTS})
  include(CTest)
  enable_testing()
  add_subdirectory(tests)
endif ()
