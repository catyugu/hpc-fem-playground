# HPC-FEM Naming Registry

This document tracks all classes, enums, structs, and namespaces in the `hpcfem` library to ensure naming uniqueness (Guideline #5).

## Namespaces

- `hpcfem` - Root namespace for all library code

## Classes

- `SolverInterface` - Abstract base class for linear solvers (src/hpcfem/solver_interface.hpp)
- `HypreAmgSolver` - Concrete HYPRE AMG solver implementation (src/hpcfem/solver_hypre_amg.hpp/cpp)

## Structs

(None yet)

## Enums

(None yet)

## Files

- `src/hpcfem/solver_interface.hpp` - Abstract solver interface
- `src/hpcfem/solver_hypre_amg.hpp` - HYPRE AMG solver header
- `src/hpcfem/solver_hypre_amg.cpp` - HYPRE AMG solver implementation

---
**Last Updated:** 2025-10-28
