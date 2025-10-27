# Benchmark targets
# These compare standard MFEM approaches vs plugin-based optimizations

# Benchmark executables
set(BENCHMARK_TARGETS
    bench_assembly
    bench_solver
    bench_combined
)

# Create each benchmark executable
foreach(target ${BENCHMARK_TARGETS})
    add_executable(${target} ${target}.cpp)
    target_link_libraries(${target} PRIVATE HPC_FEM_PLAYGROUND_LIB)
    target_include_directories(${target} PRIVATE 
        ${HPC_FEM_PLAYGROUND_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/src
    )
    
    # Set output directory
    set_target_properties(${target} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/benchmark
    )
endforeach()

# Add custom target to build all benchmarks
add_custom_target(benchmarks DEPENDS ${BENCHMARK_TARGETS})

# Installation
install(TARGETS ${BENCHMARK_TARGETS}
    RUNTIME DESTINATION bin
)

# Print benchmark info
message(STATUS "Configured benchmarks:")
foreach(target ${BENCHMARK_TARGETS})
    message(STATUS "  - ${target}")
endforeach()
