##
# cmake/ConfigureMFEM.cmake
# Encapsulate MFEM and third-party vendoring and feature toggles.
# - MFEM can be built standalone.
# - We attempt to enable MPI by default, but MPI requires OpenMPI, HYPRE and METIS.
#   If any of those are missing or fail to vendor, MPI support will be disabled.
# - SuiteSparse/UMFPACK support is only enabled when MPI dependencies are all satisfied.
##

if(NOT DEFINED VENDOR_THIRD_PARTY)
    set(VENDOR_THIRD_PARTY ON CACHE BOOL "Clone third-party deps into share" FORCE)
endif()

if(NOT DEFINED THIRD_PARTY_DIR)
    set(THIRD_PARTY_DIR ${PROJECT_SOURCE_DIR}/share CACHE PATH "Path to store third-party sources")
endif()

# MFEM feature toggles (can be overridden by environment or CI)
if(NOT DEFINED ENABLE_MFEM_MPI)
    option(ENABLE_MFEM_MPI "Enable MPI support for MFEM" ON)
else()
    set(ENABLE_MFEM_MPI ${ENABLE_MFEM_MPI})
endif()
if(NOT DEFINED ENABLE_MFEM_HYPRE)
    option(ENABLE_MFEM_HYPRE "Enable HYPRE support for MFEM" ON)
else()
    set(ENABLE_MFEM_HYPRE ${ENABLE_MFEM_HYPRE})
endif()
if(NOT DEFINED ENABLE_MFEM_SUITESPARSE)
    option(ENABLE_MFEM_SUITESPARSE "Enable SuiteSparse/UMFPACK for MFEM" ON)
else()
    set(ENABLE_MFEM_SUITESPARSE ${ENABLE_MFEM_SUITESPARSE})
endif()

# Working variables that MFEM's build will read
set(MFEM_USE_MPI ${ENABLE_MFEM_MPI})
set(MFEM_USE_HYPRE ${ENABLE_MFEM_HYPRE})
set(MFEM_USE_SUITESPARSE ${ENABLE_MFEM_SUITESPARSE})

if(VENDOR_THIRD_PARTY)
    file(MAKE_DIRECTORY ${THIRD_PARTY_DIR})
endif()

# If MPI requested, ensure we have OpenMPI, HYPRE and METIS available.
set(_mpi_deps_ok TRUE)
if(ENABLE_MFEM_MPI)
    if(NOT VENDOR_THIRD_PARTY)
        # We're not vendoring; assume system packages will be used. Keep MPI requested.
        message(STATUS "ENABLE_MFEM_MPI is ON and vendoring disabled: expecting system OpenMPI/HYPRE/METIS")
    else()
        message(STATUS "ENABLE_MFEM_MPI is ON: attempting to vendor OpenMPI, HYPRE and METIS into ${THIRD_PARTY_DIR}")

        # METIS
        set(METIS_SRC ${THIRD_PARTY_DIR}/metis)
        if(NOT EXISTS ${METIS_SRC})
            message(STATUS "Cloning METIS into ${METIS_SRC}")
            execute_process(COMMAND git clone --depth 1 https://github.com/KarypisLab/METIS.git ${METIS_SRC}
                            --branch v5.2.1    RESULT_VARIABLE _gitret_metis
                            OUTPUT_QUIET ERROR_QUIET)
            if(NOT _gitret_metis EQUAL 0)
                message(WARNING "Failed to clone METIS into ${METIS_SRC}")
            endif()
        endif()
        if(NOT EXISTS ${METIS_SRC})
            message(WARNING "METIS not available; disabling MFEM MPI support")
            set(_mpi_deps_ok FALSE)
        endif()

        message(STATUS "Looking for system-installed HYPRE")
        find_package(HYPRE QUIET)
        if(HYPRE_FOUND)
            message(STATUS "Found HYPRE via find_package")
            if(DEFINED HYPRE_INCLUDE_DIR AND NOT "${HYPRE_INCLUDE_DIR}" STREQUAL "")
                set(HYPRE_INCLUDE_DIR ${HYPRE_INCLUDE_DIR})
            elseif(DEFINED HYPRE_INCLUDE_DIRS AND NOT "${HYPRE_INCLUDE_DIRS}" STREQUAL "")
                set(HYPRE_INCLUDE_DIR ${HYPRE_INCLUDE_DIRS})
            endif()
        else()
            message(WARNING "HYPRE not found on the system; disabling MFEM MPI support")
            set(_mpi_deps_ok FALSE)
        endif()

        # OpenMPI (source clone only — building is left to the developer/CI)
        set(OPENMPI_SRC ${THIRD_PARTY_DIR}/openmpi)
        if(NOT EXISTS ${OPENMPI_SRC})
            message(STATUS "Cloning OpenMPI into ${OPENMPI_SRC}")
            execute_process(COMMAND git clone --depth 1 https://github.com/open-mpi/ompi.git ${OPENMPI_SRC}
                            --branch v4.1.1    RESULT_VARIABLE _gitret_ompi
                            OUTPUT_QUIET ERROR_QUIET)
            if(NOT _gitret_ompi EQUAL 0)
                message(WARNING "Failed to clone OpenMPI into ${OPENMPI_SRC}")
            endif()
        endif()
        if(NOT EXISTS ${OPENMPI_SRC})
            message(WARNING "OpenMPI not available; disabling MFEM MPI support")
            set(_mpi_deps_ok FALSE)
        else()
            set(OPENMPI_INCLUDE_DIR ${OPENMPI_SRC}/opal/include)
        endif()

        if(NOT _mpi_deps_ok)
            message(STATUS "MPI dependencies incomplete; turning off MFEM MPI support")
            set(MFEM_USE_MPI OFF)
            set(MFEM_USE_HYPRE OFF)
            set(MFEM_USE_SUITESPARSE OFF)
        else()
            message(STATUS "All MPI vendored dependencies present: MFEM MPI support remains enabled")
        endif()
    endif()
endif()

# SuiteSparse/UMFPACK should only be enabled when MPI deps are satisfied
if(NOT MFEM_USE_MPI)
    if(MFEM_USE_SUITESPARSE)
        message(STATUS "Disabling SuiteSparse/UMFPACK because MPI support is not available")
        set(MFEM_USE_SUITESPARSE OFF)
    endif()
endif()

if(VENDOR_THIRD_PARTY AND MFEM_USE_SUITESPARSE)
    # SuiteSparse (for UMFPACK)
    set(SUITESPARSE_SRC ${THIRD_PARTY_DIR}/SuiteSparse)
    if(NOT EXISTS ${SUITESPARSE_SRC})
        message(STATUS "Cloning SuiteSparse into ${SUITESPARSE_SRC}")
        execute_process(COMMAND git clone --depth 1 https://github.com/DrTimothyAldenDavis/SuiteSparse.git ${SUITESPARSE_SRC}
                        --branch v7.9.0  RESULT_VARIABLE _gitret_suitesparse
                        OUTPUT_QUIET ERROR_QUIET)
        if(NOT _gitret_suitesparse EQUAL 0)
            message(WARNING "Failed to clone SuiteSparse into ${SUITESPARSE_SRC}")
        endif()
    endif()
    if(EXISTS ${SUITESPARSE_SRC})
        set(UMF_INCLUDE_DIR ${SUITESPARSE_SRC}/UMFPACK/Include)
    endif()
endif()

message(STATUS "Third-party dir: ${THIRD_PARTY_DIR}")

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

# Expose variables to parent / other CMakeLists
set(MFEM_INCLUDE_DIR ${MFEM_INCLUDE_DIR} CACHE PATH "Path to MFEM sources")
set(MFEM_USE_MPI ${MFEM_USE_MPI} CACHE BOOL "Whether to build MFEM with MPI")
set(MFEM_USE_HYPRE ${MFEM_USE_HYPRE} CACHE BOOL "Whether to build MFEM with HYPRE")
set(MFEM_USE_SUITESPARSE ${MFEM_USE_SUITESPARSE} CACHE BOOL "Whether to build MFEM with SuiteSparse/UMFPACK")
