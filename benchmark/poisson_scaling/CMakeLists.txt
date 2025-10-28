# Poisson scaling benchmark

add_executable(benchmark_poisson main.cpp)

# Link against hpcfem library
target_link_libraries(benchmark_poisson PRIVATE hpcfem)

# Set C++17 standard
target_compile_features(benchmark_poisson PRIVATE cxx_std_17)

# One-level DDM scaling benchmark (Phase 2.1)
add_executable(benchmark_poisson_ddm1 poisson_scaling_ddm1.cpp)
target_link_libraries(benchmark_poisson_ddm1 PRIVATE hpcfem)
target_compile_features(benchmark_poisson_ddm1 PRIVATE cxx_std_17)

# MPI support
if(MFEM_USE_MPI)
    find_package(MPI REQUIRED)
    target_link_libraries(benchmark_poisson PRIVATE MPI::MPI_CXX)
    target_link_libraries(benchmark_poisson_ddm1 PRIVATE MPI::MPI_CXX)
endif()
