cmake_minimum_required(VERSION 3.16)
project(HPC_FEM_PLAYGROUND CXX)

# Options to vendor third-party dependencies into `share`
option(VENDOR_THIRD_PARTY "Clone third-party deps into share" ON)
set(THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/share CACHE PATH "Path to store third-party sources")

# MFEM feature toggles (can be overridden by environment or CI)
option(ENABLE_MFEM_MPI "Enable MPI support for MFEM" ON)
option(ENABLE_MFEM_HYPRE "Enable HYPRE support for MFEM" ON)
option(ENABLE_MFEM_SUITESPARSE "Enable SuiteSparse/UMFPACK for MFEM" ON)

set(MFEM_USE_MPI ${ENABLE_MFEM_MPI})
set(MFEM_USE_HYPRE ${ENABLE_MFEM_HYPRE})
set(MFEM_USE_SUITESPARSE ${ENABLE_MFEM_SUITESPARSE})

# If requested, vendor (git-clone) third-party sources into share/third_party
include(configurations/ConfigureMFEM.cmake)

add_subdirectory(${MFEM_INCLUDE_DIR} ${CMAKE_BINARY_DIR}/mfem_build)
set(HPC_FEM_PLAYGROUND_INCLUDE_DIR
    ${MFEM_INCLUDE_DIR}
)
# Make a lib that links to MFEM
add_library(HPC_FEM_PLAYGROUND_LIB INTERFACE)
target_include_directories(HPC_FEM_PLAYGROUND_LIB INTERFACE ${HPC_FEM_PLAYGROUND_INCLUDE_DIR})
target_link_libraries(HPC_FEM_PLAYGROUND_LIB INTERFACE mfem)

add_subdirectory(src)
add_subdirectory(example)
add_subdirectory(benchmark)