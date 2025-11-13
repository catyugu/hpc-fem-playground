cmake_minimum_required(VERSION 3.16)
project(HPC_FEM_PLAYGROUND VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Include CPM for dependency management
include(cmake/CPM.cmake)

# By Default, enable MFEM use of MPI
set(MFEM_USE_MPI ON CACHE BOOL "Enable MPI support in MFEM" FORCE)
if (MFEM_USE_MPI)
    message(STATUS "Configuring MFEM with MPI support")
    find_package(MPI REQUIRED)
    set(MFEM_USE_METIS_5 ON)
    set(MFEM_USE_HYPRE ON)
    set(MFEM_USE_SUITESPARSE ON)
    set(HYPRE_INCLUDE_DIRS /usr/include/hypre)
    set(HYPRE_LIBRARIES /usr/lib/x86_64-linux-gnu/libHYPRE.so)
endif()

if(NOT DEFINED THIRD_PARTY_DIR)
    set(THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/share CACHE PATH "Path to store third-party sources")
endif()

# MFEM source (vendor or clone)
set(MFEM_INCLUDE_DIR ${THIRD_PARTY_DIR}/mfem)
if(NOT EXISTS ${MFEM_INCLUDE_DIR})
    message(STATUS "MFEM not found locally, downloading from remote...")
    execute_process(
        COMMAND git clone --depth 1 --branch v4.8 https://github.com/mfem/mfem.git ${MFEM_INCLUDE_DIR}
        RESULT_VARIABLE GIT_CLONE_RESULT
        OUTPUT_QUIET ERROR_QUIET
    )
    if(NOT GIT_CLONE_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to clone MFEM repository")
    endif()
    message(STATUS "MFEM download completed")
endif()
# Prevent MFEM from configuring/building its own examples and tests
# Set these cache variables before adding MFEM so they take effect.
set(MFEM_ENABLE_EXAMPLES OFF CACHE BOOL "Disable building MFEM examples" FORCE)
set(MFEM_ENABLE_TESTING  OFF CACHE BOOL "Disable building MFEM tests" FORCE)

add_subdirectory(${MFEM_INCLUDE_DIR})

# Build the hpcfem library
add_subdirectory(src)

# Enable testing
enable_testing()
add_subdirectory(tests)

add_subdirectory(example)
add_subdirectory(benchmark)