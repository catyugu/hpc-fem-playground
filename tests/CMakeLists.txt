# Test configuration for HPC-FEM-PLAYGROUND
cmake_minimum_required(VERSION 3.16)

enable_testing()

# Test 1: Verify MFEM linkage
add_executable(test_mfem_linkage test_mfem_linkage.cpp)
target_link_libraries(test_mfem_linkage PRIVATE mfem)
add_test(NAME test_mfem_linkage COMMAND test_mfem_linkage)

# Test 2: MFEM baseline example (Phase 0, Step 0.2)
add_executable(test_mfem_baseline test_mfem_baseline.cpp)
target_link_libraries(test_mfem_baseline PRIVATE mfem)
add_test(NAME test_mfem_baseline COMMAND test_mfem_baseline)

# Test 3: Electrostatics MMS test (Phase 1, Step 1.1)
add_executable(test_physics_electrostatics test_physics_electrostatics.cpp)
target_link_libraries(test_physics_electrostatics PRIVATE mfem)
add_test(NAME test_physics_electrostatics COMMAND test_physics_electrostatics)

# Test 4: Thermal MMS test (Phase 1, Step 1.2)
add_executable(test_physics_thermal test_physics_thermal.cpp)
target_link_libraries(test_physics_thermal PRIVATE mfem)
add_test(NAME test_physics_thermal COMMAND test_physics_thermal)

# For parallel tests with MPI
if(MFEM_USE_MPI)
    find_package(MPI REQUIRED)
    add_test(NAME test_mfem_baseline_parallel 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_mfem_baseline> ${MPIEXEC_POSTFLAGS})
    add_test(NAME test_physics_electrostatics_parallel 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_physics_electrostatics> ${MPIEXEC_POSTFLAGS})
    add_test(NAME test_physics_thermal_parallel 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_physics_thermal> ${MPIEXEC_POSTFLAGS})
endif()
