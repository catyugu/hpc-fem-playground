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

# Test 5: One-way coupled Joule heating (Phase 1, Step 1.3)
add_executable(test_physics_coupling test_physics_coupling.cpp)
target_link_libraries(test_physics_coupling PRIVATE mfem)
add_test(NAME test_physics_coupling COMMAND test_physics_coupling)

# Test 6: Solver interface compile test (Phase 2, Step 2.1)
add_executable(test_solver_interface test_solver_interface.cpp)
target_link_libraries(test_solver_interface PRIVATE mfem)
add_test(NAME test_solver_interface COMMAND test_solver_interface)

# Test 7: HYPRE AMG Solver test (Phase 2, Step 2.2)
add_executable(test_solver_hypre_amg 
    test_solver_hypre_amg.cpp
    ../src/hpcfem/solver_hypre_amg.cpp)
target_include_directories(test_solver_hypre_amg PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_solver_hypre_amg PRIVATE mfem)
add_test(NAME test_solver_hypre_amg COMMAND test_solver_hypre_amg)

# Test 8: Problem abstraction test (Phase 2, Step 2.3)
add_executable(test_problem_abstraction 
    test_problem_abstraction.cpp
    ../src/hpcfem/solver_hypre_amg.cpp
    ../src/hpcfem/physics_electrostatics.cpp
    ../src/hpcfem/fem_problem.cpp)
target_include_directories(test_problem_abstraction PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_link_libraries(test_problem_abstraction PRIVATE mfem)
add_test(NAME test_problem_abstraction COMMAND test_problem_abstraction)

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
    add_test(NAME test_physics_coupling_parallel 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_physics_coupling> ${MPIEXEC_POSTFLAGS})
    add_test(NAME test_solver_interface_parallel 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_solver_interface> ${MPIEXEC_POSTFLAGS})
    add_test(NAME test_solver_hypre_amg_parallel 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_solver_hypre_amg> ${MPIEXEC_POSTFLAGS})
    add_test(NAME test_problem_abstraction_parallel 
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 
             ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test_problem_abstraction> ${MPIEXEC_POSTFLAGS})
endif()
