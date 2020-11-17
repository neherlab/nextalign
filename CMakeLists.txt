cmake_minimum_required(VERSION 3.14)
cmake_policy(SET CMP0069 NEW)

include(packages/cmake/quiet.cmake)
include(packages/cmake/CoolCFlags.cmake)
include(packages/cmake/CoolCXXFlags.cmake)

project(nextalign_all)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_INSTALL_MESSAGE LAZY)

add_subdirectory(packages/nextalign)
add_subdirectory(packages/nextalign_cli)
