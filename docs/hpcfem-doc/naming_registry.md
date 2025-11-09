# HPC-FEM Naming Registry

This document tracks all classes, enums, structs, and namespaces in the `hpcfem` library to ensure naming uniqueness (Guideline #5).

**Last Updated:** 2025-10-28

## Namespaces

- `hpcfem` - Root namespace for all library code

## Core Library Classes

### Abstract Interfaces

- `PhysicsInterface` - Abstract base class for physics modules (src/hpcfem/physics_interface.hpp)
  - **Status:** Pure virtual ✓
  - **Purpose:** Defines contract for physics assembly
  - **Key Methods:** `assemble()`, `getFiniteElementSpace()`
  
- `SolverInterface` - Abstract base class for linear solvers (src/hpcfem/solver_interface.hpp)
  - **Status:** Pure virtual ✓
  - **Purpose:** Defines contract for linear system solvers
  - **Key Methods:** `solve()`

### Concrete Implementations

- `PhysicsCavityEigen` - 3D Maxwell eigenvalue problem solver (src/hpcfem/physics/physics_cavity_eigen.hpp/cpp)
  - **Purpose:** Solves ∇×∇×E = λE for electromagnetic cavity modes
  - **Equation:** Maxwell eigenvalue problem with PEC boundaries
  - **Key Feature:** Uses Nedelec elements and HypreAMS preconditioner
  - **Application:** Waveguide mode analysis, S-parameter port definitions
  
- `HypreAmgSolver` - HYPRE AMG solver wrapper (src/hpcfem/solver_hypre_amg.hpp/cpp)
  - **Implements:** `SolverInterface`
  - **Backend:** HYPRE BoomerAMG

## Test Classes (Not in Library)

### Mock Implementations (tests/)

- `MockSolver` - Mock implementation of `SolverInterface` for testing (test_solver_interface_mock.cpp)
- `MockPhysics` - Mock implementation of `PhysicsInterface` for testing (test_physics_interface_mock.cpp)

## Files by Category

### Core Interfaces

- `SolverInterface` - Abstract base class for linear solvers (src/hpcfem/core/solver_interface.hpp)

- `src/hpcfem/solver_interface.hpp` - Abstract solver interface

### Physics Modules

- `src/hpcfem/physics/physics_cavity_eigen.hpp` - 3D Maxwell eigenvalue solver header
- `src/hpcfem/physics/physics_cavity_eigen.cpp` - 3D Maxwell eigenvalue solver implementation
- `src/hpcfem/physics/physics_waveguide_eigen.hpp` - 2D waveguide eigenvalue solver header
- `src/hpcfem/physics/physics_waveguide_eigen.cpp` - 2D waveguide eigenvalue solver implementation

### Solver Modules

- `src/hpcfem/core/solver_interface.hpp` - Abstract solver interface
- `src/hpcfem/solvers/solver_hypre_amg.hpp` - HYPRE AMG solver header
- `src/hpcfem/solvers/solver_hypre_amg.cpp` - HYPRE AMG solver implementation

### Tests

- `tests/test_solver_interface.cpp` - Solver interface compilation test
- `tests/test_solver_interface_mock.cpp` - Solver interface contract test
- `tests/test_physics_interface_mock.cpp` - Physics interface contract test
- `tests/test_solver_hypre_amg.cpp` - HYPRE AMG solver functional test
- `tests/test_physics_cavity_eigen.cpp` - 3D Maxwell eigenvalue problem test
- `tests/test_physics_waveguide_eigen.cpp` - 2D waveguide eigenvalue problem test

### Examples

- `example/ex0_cylinder_dirichlet.cpp` - Simple cylinder electrostatics
- `example/ex1_cube_mixed_bc.cpp` - Cube with mixed boundary conditions
- `example/ex2_rect_waveguide.cpp` - Rectangular waveguide
- `example/ex2p_rect_waveguide_hypre.cpp` - Rectangular waveguide with HYPRE
- `example/ex3_transient_heat.cpp` - Transient heat diffusion
- `example/ex4_eigenvalue_laplacian.cpp` - Eigenvalue problem
- `example/ex5_coupled_problem.cpp` - Coupled electro-thermal problem
- `example/ex6_multiple_materials.cpp` - Nonlinear Joule heating with Q = σ|∇V|²

