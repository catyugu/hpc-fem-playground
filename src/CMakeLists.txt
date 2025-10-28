# hpcfem library
add_library(hpcfem
    hpcfem/fem_problem.cpp
    hpcfem/physics_electrostatics.cpp
    hpcfem/solver_hypre_amg.cpp
    hpcfem/solver_one_level_schwarz.cpp
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