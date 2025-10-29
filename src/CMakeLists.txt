# hpcfem library
add_library(hpcfem
    hpcfem/fem_problem.cpp
    hpcfem/physics/physics_electrostatics.cpp
    hpcfem/physics/physics_thermal.cpp
    hpcfem/physics/physics_joule_heating.cpp
    hpcfem/solvers/solver_hypre_amg.cpp
    hpcfem/solvers/solver_block_gauss_seidel.cpp
)

# Set include directories - src is the root
target_include_directories(hpcfem PUBLIC 
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link dependencies
target_link_libraries(hpcfem PUBLIC 
    mfem
)

# Export compile features
target_compile_features(hpcfem PUBLIC cxx_std_17)