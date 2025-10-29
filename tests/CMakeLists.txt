# Test configuration for HPC-FEM-PLAYGROUND
cmake_minimum_required(VERSION 3.16)

# Fetch GoogleTest using CPM
CPMAddPackage(
    NAME googletest
    GITHUB_REPOSITORY google/googletest
    VERSION 1.14.0
    OPTIONS
        "INSTALL_GTEST OFF"
        "gtest_force_shared_crt ON"
)

enable_testing()

# Automatically collect all test_*.cpp files
file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/test_*.cpp")

# Create test executables and register them with CTest
foreach(TEST_SOURCE ${TEST_SOURCES})
    # Extract test name from file path (e.g., test_solver_interface.cpp -> test_solver_interface)
    get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
    
    # Create executable
    add_executable(${TEST_NAME} ${TEST_SOURCE})
    
    # Link libraries (all tests use hpcfem and gtest_main)
    target_link_libraries(${TEST_NAME} PRIVATE hpcfem gtest_main)
    
    # Register with CTest for serial execution
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
    
    # Add parallel MPI variant if MPI is available
    if(MFEM_USE_MPI)
        add_test(NAME ${TEST_NAME}_parallel 
                 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
                 ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${TEST_NAME}> ${MPIEXEC_POSTFLAGS})
    endif()
endforeach()

# Find MPI package if needed (referenced above)
if(MFEM_USE_MPI)
    find_package(MPI REQUIRED)
endif()
