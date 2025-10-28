# HPC-FEM Naming Registry

This document tracks all classes, enums, structs, and namespaces in the `hpcfem` library to ensure naming uniqueness (Guideline #5).

**Last Updated:** 2025-10-28 (Phase 3.1 - In Progress)

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
  - **Phase:** 3.1
  
- `HypreAmgSolver` - HYPRE AMG solver wrapper (src/hpcfem/solver_hypre_amg.hpp/cpp)
  - **Implements:** `SolverInterface`
  - **Backend:** HYPRE BoomerAMG

- `OneLevelSchwarz` - One-level overlapping Schwarz preconditioner (src/hpcfem/solver_one_level_schwarz.hpp/cpp)
  - **Extends:** `mfem::Solver`
  - **Purpose:** Domain decomposition preconditioner (additive Schwarz)

- `TwoLevelSchwarz` - Two-level overlapping Schwarz preconditioner with AMG coarse solver (src/hpcfem/solver_two_level_schwarz.hpp/cpp)
  - **Extends:** `mfem::Solver`
  - **Purpose:** Framework for scalable two-level DD method (coarse correction TODO)

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
- `src/hpcfem/physics_interface.hpp` - Abstract physics interface
- `src/hpcfem/solver_interface.hpp` - Abstract solver interface

### Physics Modules
- `src/hpcfem/physics_electrostatics.hpp` - Electrostatics physics header
- `src/hpcfem/physics_electrostatics.cpp` - Electrostatics physics implementation
- `src/hpcfem/physics_thermal.hpp` - Thermal diffusion physics header (Phase 3.1)
- `src/hpcfem/physics_thermal.cpp` - Thermal diffusion physics implementation (Phase 3.1)

### Solver Modules
- `src/hpcfem/solver_hypre_amg.hpp` - HYPRE AMG solver header
- `src/hpcfem/solver_hypre_amg.cpp` - HYPRE AMG solver implementation
- `src/hpcfem/solver_one_level_schwarz.hpp` - One-level Schwarz preconditioner header (Phase 2.1)
- `src/hpcfem/solver_one_level_schwarz.cpp` - One-level Schwarz preconditioner implementation (Phase 2.1)
- `src/hpcfem/solver_two_level_schwarz.hpp` - Two-level Schwarz preconditioner header (Phase 2.2)
- `src/hpcfem/solver_two_level_schwarz.cpp` - Two-level Schwarz preconditioner implementation (Phase 2.2)

### Problem Orchestration
- `src/hpcfem/fem_problem.hpp` - FEM problem orchestrator header
- `src/hpcfem/fem_problem.cpp` - FEM problem orchestrator implementation

### Tests
- `tests/test_solver_interface.cpp` - Solver interface compilation test
- `tests/test_solver_interface_mock.cpp` - Solver interface contract test (Phase 1.2)
- `tests/test_physics_interface_mock.cpp` - Physics interface contract test (Phase 1.2)
- `tests/test_solver_hypre_amg.cpp` - HYPRE AMG solver functional test
- `tests/test_problem_abstraction.cpp` - End-to-end abstraction test
- `tests/test_physics_electrostatics.cpp` - Electrostatics MMS test
- `tests/test_physics_thermal.cpp` - Thermal physics test
- `tests/test_physics_coupling.cpp` - Coupled physics test
- `tests/test_ddm_one_level.cpp` - One-level Schwarz correctness test (Phase 2.1)
- `tests/test_ddm_two_level.cpp` - Two-level Schwarz correctness test (Phase 2.2)

### Benchmarks
- `benchmark/poisson_scaling/main.cpp` - Original Poisson scaling benchmark
- `benchmark/poisson_scaling/poisson_scaling_ddm1.cpp` - One-level Schwarz scaling benchmark (Phase 2.1)
- `benchmark/poisson_scaling/poisson_scaling_ddm2.cpp` - Two-level Schwarz framework benchmark (Phase 2.2)

---
