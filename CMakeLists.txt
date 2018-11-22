CMAKE_MINIMUM_REQUIRED(VERSION 3.9)

# required for unit testing using CMake's ctest command
include(CTest)
enable_testing()

project(lifev-pyorb VERSION 1.0.0 DESCRIPTION "lifev-pyorb description")

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)
include(ExternalProject)


set(CMAKE_BUILD_TYPE Release)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include(TPLs/lifev.cmake)

add_subdirectory(include)
add_subdirectory(TPLs)
add_subdirectory(sources)
add_subdirectory(test)

add_library(lifev-pyorb SHARED sources/test.cpp )


