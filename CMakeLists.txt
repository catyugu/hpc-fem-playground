cmake_minimum_required(VERSION 3.16)
project(HPC_FEM_PLAYGROUND CXX)

set(MFEM_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/share/mfem)
if(NOT EXISTS ${MFEM_INCLUDE_DIR})
    message(STATUS "MFEM not found locally, downloading from remote...")
    execute_process(
        COMMAND git clone --depth 1 --branch v4.8 https://github.com/mfem/mfem.git ${MFEM_INCLUDE_DIR}
        RESULT_VARIABLE GIT_CLONE_RESULT
    )
    if(NOT GIT_CLONE_RESULT EQUAL 0)
        message(FATAL_ERROR "Failed to clone MFEM repository")
    endif()
    message(STATUS "MFEM download completed")
endif()

# If MFEM was cloned into share/mfem, build it as a subdirectory so the
# "mfem" CMake target is available for linking. We build it in a separate
# build subdirectory under the main build dir to avoid in-source MFEM builds.
if (EXISTS ${MFEM_INCLUDE_DIR}/CMakeLists.txt)
    message(STATUS "Adding MFEM subdirectory: ${MFEM_INCLUDE_DIR}")
    add_subdirectory(${MFEM_INCLUDE_DIR} ${CMAKE_BINARY_DIR}/mfem_build)
else()
    message(WARNING "MFEM CMakeLists.txt not found in ${MFEM_INCLUDE_DIR}; expecting system-installed libmfem")
endif()

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