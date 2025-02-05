cmake_minimum_required(VERSION 3.24.0)

project(ale CXX)

find_package(MPI REQUIRED)

# Require openMP for this project
include(FindOpenMP)
if(OPENMP_FOUND)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else(OPENMP_FOUND)
  message("ERROR: OpenMP could not be found on your system.")
endif(OPENMP_FOUND)

set(BOOST_INCLUDE_LIBRARIES algorithm math mpi serialization timer)
set(BOOST_ENABLE_MPI ON)

# Fetching dependencies can be quite long, so we should log the progress
set(FETCHCONTENT_QUIET false)
include(FetchContent)

# Downloading as tar.gz archive is way faster than through git
FetchContent_Declare(
    Boost
    URL https://github.com/boostorg/boost/releases/download/boost-1.86.0/boost-1.86.0-cmake.tar.gz
    URL_HASH SHA256=c62ce6e64d34414864fef946363db91cea89c1b90360eabed0515f0eda74c75c
    FIND_PACKAGE_ARGS COMPONENTS "${BOOST_INCLUDE_LIBRARIES}" CONFIG
)

FetchContent_MakeAvailable(Boost)
set(PATCH_DIR "${PROJECT_SOURCE_DIR}/patches")

set(BPP_CORE_PATCH_COMMAND git apply "${PATCH_DIR}/bpp-core.patch" || true)
set(BPP_SEQ_PATCH_COMMAND git apply "${PATCH_DIR}/bpp-seq.patch" || true)
set(BPP_PHYL_PATCH_COMMAND git apply "${PATCH_DIR}/bpp-phyl.patch" || true)

FetchContent_Declare(
    bpp-core
    GIT_REPOSITORY https://github.com/BioPP/bpp-core
    GIT_TAG v2.4.1
    PATCH_COMMAND ${BPP_CORE_PATCH_COMMAND}
    OVERRIDE_FIND_PACKAGE
)

FetchContent_Declare(
    bpp-seq
    GIT_REPOSITORY https://github.com/BioPP/bpp-seq
    GIT_TAG v2.4.1
    PATCH_COMMAND ${BPP_SEQ_PATCH_COMMAND}
    OVERRIDE_FIND_PACKAGE
)

FetchContent_Declare(
    bpp-phyl
    GIT_REPOSITORY https://github.com/BioPP/bpp-phyl
    GIT_TAG v2.4.1
    PATCH_COMMAND ${BPP_PHYL_PATCH_COMMAND}
)

FetchContent_MakeAvailable(bpp-core bpp-seq bpp-phyl)

set(BPP_LIBRARIES bpp-core-static bpp-seq-static bpp-phyl-static)

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/lib")

add_subdirectory(src)