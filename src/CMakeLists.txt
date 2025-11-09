# hpcfem library
add_library(hpcfem
    hpcfem/physics/physics_cavity_eigen.cpp
    hpcfem/physics/physics_waveguide_eigen.cpp
    hpcfem/solvers/solver_hypre_amg.cpp
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