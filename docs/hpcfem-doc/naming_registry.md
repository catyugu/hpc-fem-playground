# HPC-FEM Naming Registry

This document tracks all classes, enums, structs, and namespaces in the `hpcfem` library to ensure naming uniqueness (Guideline #5).

## Namespaces

- `hpcfem` - Root namespace for all library code

## Classes

- `SolverInterface` - Abstract base class for linear solvers (src/hpcfem/solver_interface.hpp)
- `HypreAmgSolver` - Concrete HYPRE AMG solver implementation (src/hpcfem/solver_hypre_amg.hpp/cpp)
- `DdmSchwarzSolver` - Domain decomposition Schwarz solver (src/hpcfem/solver_ddm_schwarz.hpp/cpp)
- `PhysicsInterface` - Abstract base class for physics modules (src/hpcfem/physics_interface.hpp)
- `ElectrostaticsPhysics` - Concrete electrostatics physics implementation (src/hpcfem/physics_electrostatics.hpp/cpp)
- `FemProblem` - Top-level FEM problem orchestrator (src/hpcfem/fem_problem.hpp/cpp)

## Structs

(None yet)

## Enums

(None yet)

## Files

- `src/hpcfem/solver_interface.hpp` - Abstract solver interface
- `src/hpcfem/solver_hypre_amg.hpp` - HYPRE AMG solver header
- `src/hpcfem/solver_hypre_amg.cpp` - HYPRE AMG solver implementation
- `src/hpcfem/solver_ddm_schwarz.hpp` - DDM Schwarz solver header
- `src/hpcfem/solver_ddm_schwarz.cpp` - DDM Schwarz solver implementation
- `src/hpcfem/physics_interface.hpp` - Abstract physics interface
- `src/hpcfem/physics_electrostatics.hpp` - Electrostatics physics header
- `src/hpcfem/physics_electrostatics.cpp` - Electrostatics physics implementation
- `src/hpcfem/fem_problem.hpp` - FEM problem orchestrator header
- `src/hpcfem/fem_problem.cpp` - FEM problem orchestrator implementation

---
**Last Updated:** 2025-10-28
