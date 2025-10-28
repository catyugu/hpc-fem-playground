# Joule Heating Solvers Benchmark
add_executable(benchmark_joule_heating_solvers main.cpp)

target_link_libraries(benchmark_joule_heating_solvers
    PRIVATE
    hpcfem
    mfem
)

# Enable MPI if available
if(MFEM_USE_MPI)
    target_link_libraries(benchmark_joule_heating_solvers PRIVATE MPI::MPI_CXX)
endif()
