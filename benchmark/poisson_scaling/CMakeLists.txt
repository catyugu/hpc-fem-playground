# Poisson scaling benchmark

# Source files from hpcfem
set(HPCFEM_SOURCES
    ${CMAKE_SOURCE_DIR}/src/hpcfem/solver_hypre_amg.cpp
    ${CMAKE_SOURCE_DIR}/src/hpcfem/physics_electrostatics.cpp
    ${CMAKE_SOURCE_DIR}/src/hpcfem/fem_problem.cpp
)

add_executable(benchmark_poisson main.cpp ${HPCFEM_SOURCES})

# Link against MFEM
target_link_libraries(benchmark_poisson PRIVATE mfem)

# Set C++17 standard
target_compile_features(benchmark_poisson PRIVATE cxx_std_17)

# Include directories - no relative paths
target_include_directories(benchmark_poisson PRIVATE 
    ${CMAKE_SOURCE_DIR}/src
)

# MPI support
if(MFEM_USE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(benchmark_poisson PRIVATE MPI::MPI_CXX)
endif()

