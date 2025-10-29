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

- `ElectrostaticsPhysics` - Electrostatics physics (Poisson equation) (src/hpcfem/physics_electrostatics.hpp/cpp)
  - **Implements:** `PhysicsInterface`
  - **Equation:** ∇·(σ∇φ) = f

- `ThermalPhysics` - Thermal diffusion physics (src/hpcfem/physics_thermal.hpp/cpp)
  - **Implements:** `PhysicsInterface`
  - **Equation:** ∇·(κ∇T) = Q

- `JouleHeatingPhysics` - Coupled electro-thermal Joule heating physics (src/hpcfem/physics_joule_heating.hpp/cpp)
  - **Does NOT implement:** `PhysicsInterface` (uses BlockOperator, specialized interface)
  - **Equation:** 2×2 block system `[K_e 0; C K_t][V; T] = [f_e; f_t]`
  - **Key Feature:** Multiphysics with nonlinear Joule heating Q = σ|∇V|²
  - **Nonlinearity:** Solved via Picard iteration with voltage-dependent thermal source
  
- `JouleHeatingCoefficient` - Joule heating source term coefficient (src/hpcfem/joule_heating_coefficient.hpp)
  - **Extends:** `mfem::Coefficient`
  - **Formula:** Q = σ|∇V|²
  - **Purpose:** Evaluates nonlinear Joule heating power density at quadrature points
  
- `HypreAmgSolver` - HYPRE AMG solver wrapper (src/hpcfem/solver_hypre_amg.hpp/cpp)
  - **Implements:** `SolverInterface`
  - **Backend:** HYPRE BoomerAMG

### Orchestration

- `FemProblem` - Top-level FEM problem orchestrator (src/hpcfem/fem_problem.hpp/cpp)
  - **Purpose:** Coordinates mesh, physics, and solver
  - **Key Methods:** `assemble()`, `solve()`, `getSolution()`

## Test Classes (Not in Library)

### Mock Implementations (tests/)

- `MockSolver` - Mock implementation of `SolverInterface` for testing (test_solver_interface_mock.cpp)
- `MockPhysics` - Mock implementation of `PhysicsInterface` for testing (test_physics_interface_mock.cpp)

## Structs

(None yet)

## Enums

(None yet)

## Files by Category

### Core Interfaces

- `SolverInterface` - Abstract base class for linear solvers (src/hpcfem/core/solver_interface.hpp)

- `src/hpcfem/solver_interface.hpp` - Abstract solver interface

### Physics Modules

- `src/hpcfem/physics/physics_electrostatics.hpp` - Electrostatics physics header
- `src/hpcfem/physics/physics_electrostatics.cpp` - Electrostatics physics implementation
- `src/hpcfem/physics/physics_thermal.hpp` - Thermal diffusion physics header
- `src/hpcfem/physics/physics_thermal.cpp` - Thermal diffusion physics implementation
- `src/hpcfem/physics/physics_joule_heating.hpp` - Joule heating coupled physics header

### Solver Modules

- `src/hpcfem/core/solver_interface.hpp` - Abstract solver interface
- `src/hpcfem/solvers/solver_hypre_amg.hpp` - HYPRE AMG solver header
- `src/hpcfem/solvers/solver_hypre_amg.cpp` - HYPRE AMG solver implementation

### Helper/Coefficient Classes

- `src/hpcfem/physics/joule_heating_coefficient.hpp` - Joule heating coefficient Q = σ|∇V|²

### Problem Orchestration

- `src/hpcfem/core/fem_problem.hpp` - FEM problem orchestrator header
- `src/hpcfem/core/fem_problem.cpp` - FEM problem orchestrator implementation

### Tests

- `tests/test_solver_interface.cpp` - Solver interface compilation test
- `tests/test_solver_interface_mock.cpp` - Solver interface contract test
- `tests/test_physics_interface_mock.cpp` - Physics interface contract test
- `tests/test_solver_hypre_amg.cpp` - HYPRE AMG solver functional test
- `tests/test_problem_abstraction.cpp` - End-to-end abstraction test
- `tests/test_physics_electrostatics.cpp` - Electrostatics MMS test
- `tests/test_physics_coupling.cpp` - Coupled physics test
- `tests/test_physics_joule_heating_block.cpp` - Joule heating BlockOperator assembly test

### Examples

- `example/ex0_cylinder_dirichlet.cpp` - Simple cylinder electrostatics
- `example/ex1_cube_mixed_bc.cpp` - Cube with mixed boundary conditions
- `example/ex2_rect_waveguide.cpp` - Rectangular waveguide
- `example/ex3_transient_heat.cpp` - Transient heat diffusion
- `example/ex4_eigenvalue_laplacian.cpp` - Eigenvalue problem
- `example/ex5_coupled_problem.cpp` - Coupled electro-thermal problem
- `example/ex6_nonlinear_joule_heating.cpp` - Nonlinear Joule heating with Q = σ|∇V|²

### Benchmarks

- `benchmark/poisson_scaling/main.cpp` - Poisson scaling benchmark
- `benchmark/joule_heating_solvers/main.cpp` - Joule heating solver comparison

---
