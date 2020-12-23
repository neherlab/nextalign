cmake_minimum_required(VERSION 3.10)
cmake_policy(SET CMP0069 NEW)

project(nextalign_cli)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_MESSAGE LAZY)

find_package(boost_base 1.68.0 REQUIRED)
find_package(cxxopts 2.2.1 REQUIRED)
find_package(fmt 7.1.0 REQUIRED)
find_package(TBB REQUIRED)

if (${NEXTALIGN_STATIC_BUILD})
  set(STATIC_BUILD_FLAGS "-static -static-libstdc++ -static-libgcc")
endif ()

add_executable(${PROJECT_NAME}
  src/cli.cpp
  )

set_property(TARGET ${PROJECT_NAME} PROPERTY C_STANDARD 11)
set_property(TARGET ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

target_compile_definitions(${PROJECT_NAME}
  PRIVATE -DPROJECT_NAME="${PROJECT_NAME}"
  )

target_include_directories(${PROJECT_NAME}
  PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src"
  )

target_link_libraries(${PROJECT_NAME}
  PRIVATE
  boost_base::boost_base
  cxxopts::cxxopts
  fmt::fmt
  nextalign
  TBB::tbb
  ${STATIC_BUILD_FLAGS}
  )
