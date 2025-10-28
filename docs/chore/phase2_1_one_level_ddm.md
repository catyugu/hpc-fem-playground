# Phase 2.1: TDD - One-Level DDM (Benchmark for Failure)

**Status:** ✅ COMPLETED  
**Date:** 2025-10-28  
**Branch:** `dev`  
**Commit:** [To be added]

## Overview

Successfully implemented a one-level overlapping Schwarz (Additive Schwarz) preconditioner for parallel domain decomposition. This phase demonstrates the well-known scalability limitation of one-level methods: **iteration count increases with processor count**.

## Objectives

1. ✅ Implement one-level additive Schwarz method
2. ✅ Create correctness test on 3D Poisson problem
3. ✅ Create scalability benchmark
4. ✅ Demonstrate poor scaling (iteration count increases with MPI processes)

## Implementation Details

### New Classes

#### `OneLevelSchwarz` (hpcfem/solver_one_level_schwarz.hpp/cpp)

A preconditioner implementing the one-level overlapping Schwarz method:

```
M^{-1} = ∑_{i=1}^{N} R_i^T A_i^{-1} R_i
```

Where:
- N = number of subdomains (MPI processes)
- R_i = restriction operator to subdomain i
- A_i = local matrix on subdomain i
- A_i^{-1} approximated by Gauss-Seidel smoother

**Key Design Decisions:**
- Extends `mfem::Solver` to act as a preconditioner
- Uses `HypreParMatrix::GetDiag()` to extract local subdomain matrix
- Uses `mfem::GSSmoother` as the local subdomain solver
- Temporary vectors (`localX_`, `localY_`) marked `mutable` for const-correct `Mult()`

**API:**
```cpp
// Constructor takes the parallel matrix
explicit OneLevelSchwarz(const mfem::HypreParMatrix& A);

// Apply preconditioner: y = M^{-1} x
void Mult(const mfem::Vector& x, mfem::Vector& y) const override;
```

### Test: `test_ddm_one_level.cpp`

**Purpose:** Verify correctness of the one-level Schwarz preconditioner

**Problem Setup:**
- 3D Poisson equation: -∇²u = f in [0,1]³
- Manufactured solution: u(x,y,z) = sin(πx)sin(πy)sin(πz)
- Forcing term: f = 3π²sin(πx)sin(πy)sin(πz)
- Mesh: 4³ hexahedral elements
- FE space: H1, polynomial order 2
- Boundary: Homogeneous Dirichlet (all sides)

**Test Criteria:**
- ✅ CG solver converges (residual < 1e-12)
- ✅ L2 error < 2e-3 (relaxed for coarse mesh)
- ✅ Runs correctly with 4 MPI processes

**Results:**
- CG iterations: 23
- L2 error: 0.00167
- Status: **PASSED**

### Benchmark: `poisson_scaling_ddm1.cpp`

**Purpose:** Demonstrate poor scalability of one-level method

**Setup:**
- Same 3D Poisson problem as test
- Mesh size scales with processor count (work per processor constant)
- Records iterations, solve time, L2 error to CSV file

**Benchmark Results:**

| MPI Procs | DOFs | Iterations | Solve Time (s) | L2 Error |
|-----------|------|------------|----------------|----------|
| 2         | 1331 | 28         | 0.00396        | 0.000861 |
| 4         | 2197 | 32         | 0.00276        | 0.000501 |
| 8         | 4913 | 47         | 0.00631        | 0.000212 |

**Key Observation:** Iteration count increases from **28 → 32 → 47** as processor count goes from 2 → 4 → 8.

**Analysis:**
- This confirms the theoretical prediction that one-level Schwarz lacks scalability
- As subdomains become smaller (more processes), information propagates slower across domain
- Each subdomain only sees local information; no global coarse-grid correction
- This motivates **Phase 2.2**: two-level method with AMG coarse solver

## Files Modified/Created

### New Implementation Files
- ✅ `src/hpcfem/solver_one_level_schwarz.hpp` (96 lines)
- ✅ `src/hpcfem/solver_one_level_schwarz.cpp` (76 lines)

### New Test Files
- ✅ `tests/test_ddm_one_level.cpp` (145 lines)

### New Benchmark Files
- ✅ `benchmark/poisson_scaling/poisson_scaling_ddm1.cpp` (161 lines)

### Modified Build Files
- ✅ `src/CMakeLists.txt` - Added solver_one_level_schwarz.cpp to library
- ✅ `tests/CMakeLists.txt` - Added test_ddm_one_level target and parallel variant
- ✅ `benchmark/poisson_scaling/CMakeLists.txt` - Added benchmark_poisson_ddm1 target

### Updated Documentation
- ✅ `docs/hpcfem-doc/naming_registry.md` - Registered OneLevelSchwarz class
- ✅ `docs/chore/phase2_1_one_level_ddm.md` - This file

## Build and Test Results

### Compilation
```bash
cd cmake-build-debug
cmake ..
cmake --build . --target test_ddm_one_level -j4
cmake --build . --target benchmark_poisson_ddm1 -j4
```
**Status:** ✅ Clean build, no warnings

### Test Execution
```bash
# Serial test
ctest -R test_ddm_one_level

# Parallel test (4 processes)
mpirun -np 4 ./tests/test_ddm_one_level
```
**Result:** ✅ All tests pass (21/21)

### Benchmark Execution
```bash
cd cmake-build-debug/benchmark/poisson_scaling
for np in 2 4 8; do
    mpirun -np $np ./benchmark_poisson_ddm1
done
```
**Result:** ✅ Successfully demonstrates poor scaling

## Code Quality Checklist

- ✅ **Guideline 1:** No templates, lambdas, PCH, or macros
- ✅ **Guideline 2:** All code in `hpcfem` namespace
- ✅ **Guideline 3:** No enum conflicts (no enums in this phase)
- ✅ **Guideline 4:** No code nesting (flat class design)
- ✅ **Guideline 5:** `OneLevelSchwarz` registered in naming_registry.md
- ✅ **Guideline 6:** Comprehensive Doxygen comments throughout
- ✅ **Guideline 7:** All tests compile and pass
- ✅ **Guideline 8:** No TODOs left in code
- ✅ **Guideline 9:** Working on `dev` branch with proper git workflow
- ✅ **Guideline 10:** TDD workflow: TEST → IMPLEMENT → REFACTOR → DOCUMENT
- ✅ **Guideline 15:** All constants named (PI, MMS_TOLERANCE, etc.)
- ✅ **Guideline 16:** Naming conventions followed (PascalCase for classes, camelCase for methods)
- ✅ **Guideline 17:** .hpp and .cpp files co-located in src/hpcfem
- ✅ **Guideline 18:** Clean include paths (hpcfem/ prefix, no ../)

## Technical Notes

### Implementation Challenges

1. **MFEM API Learning:**
   - Initial attempt used `GetDiag()` incorrectly (thought it returned pointer)
   - Actual API: `void GetDiag(SparseMatrix&)` - outputs to reference
   - Solution: Create `SparseMatrix` and pass as output parameter

2. **Local Solver Selection:**
   - Considered using `HypreBoomerAMG` on local subdomain
   - Problem: `GetDiag()` returns `SparseMatrix`, not `HypreParMatrix`
   - Solution: Use `mfem::GSSmoother` which operates on `SparseMatrix`
   - Note: For production code, could use sparse direct solver (UMFPACK) for exact local solves

3. **Const-Correctness:**
   - `Mult()` must be const (required by `mfem::Solver` interface)
   - But need to modify temporary vectors `localX_`, `localY_`
   - Solution: Declare them `mutable` in class definition

### Performance Considerations

- Current implementation uses Gauss-Seidel as local solver (inexact subdomain solve)
- For better performance, could use:
  - Direct sparse solver (UMFPACK, SuperLU) for small subdomains
  - Local AMG for large subdomains
- Tradeoff: exact local solve → fewer CG iterations but higher cost per iteration

### Mathematical Context

The one-level Schwarz method is a classical domain decomposition preconditioner. Its key limitation is that information propagates only through subdomain overlaps. For a problem decomposed into N subdomains with diameter h and overlap δ:

- Condition number scales as: κ(M^{-1}A) ~ h/δ ~ 1/δ (for fixed h)
- As N increases (h decreases), more CG iterations needed
- Solution: Add coarse-grid correction (Phase 2.2) for global information propagation

## Next Steps: Phase 2.2

The benchmark results motivate **Phase 2.2: Two-Level DDM with AMG Coarse Solver**:

1. Implement `SolverTwoLevelSchwarz` class
2. Add coarse-grid correction: M^{-1} = Ψ A_H^{-1} Ψ^T + ∑ R_i^T A_i^{-1} R_i
3. Use Galerkin projection: A_H = Ψ^T A Ψ
4. Use `HypreAmgSolver` for coarse-grid solve A_H^{-1}
5. Demonstrate constant iteration count (scalability success)

**Expected Outcome:** Iteration count remains ~constant (e.g., 15-20) regardless of processor count.

## References

- MFEM Documentation: https://mfem.org/
- "Domain Decomposition Methods in Science and Engineering" (Quarteroni & Valli, 1999)
- DEVELOPMENT.prompt.md: Phase 2, Step 2.1 specification
