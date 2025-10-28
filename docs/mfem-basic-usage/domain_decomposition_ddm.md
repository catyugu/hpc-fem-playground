# Domain Decomposition Methods (DDM) in MFEM

**Author:** HPC-FEM-Playground Team  
**Date:** October 2025  
**Phase:** 2 (Scalable Single-Physics Solvers)

## Overview

Domain Decomposition Methods (DDM) are parallel solvers that partition a global problem into overlapping or non-overlapping subdomains. This guide documents lessons learned implementing one-level and two-level Schwarz methods using MFEM.

## Key MFEM Classes for DDM

### Parallel Mesh Partitioning

```cpp
#include "mfem.hpp"

// Create or load serial mesh
mfem::Mesh* serial_mesh = new mfem::Mesh("mesh.mesh");

// Partition for MPI processes
mfem::ParMesh* pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *serial_mesh);
delete serial_mesh;
```

**What happens:**
- MFEM automatically uses METIS to partition the mesh graph
- Each MPI rank gets a local subdomain with ghost/shared elements
- `pmesh->GetSharedFaceDofs()` identifies interface DOFs

### Parallel Finite Element Spaces

```cpp
// Define FE space on parallel mesh
mfem::H1_FECollection fec(order, dim);
mfem::ParFiniteElementSpace fespace(pmesh, &fec);

// Get dimensions
HYPRE_BigInt global_dofs = fespace.GlobalTrueVSize();  // Total DOFs across all ranks
int local_dofs = fespace.GetTrueVSize();                // DOFs owned by this rank
```

**True DOFs vs. Local DOFs:**
- **Local DOFs**: All DOFs visible to this rank (owned + shared)
- **True DOFs**: DOFs uniquely owned by this rank (no duplicates)
- MFEM handles the parallel communication automatically

### Extracting Local Subdomain Matrices

For implementing Schwarz preconditioners, you need the local subdomain operator:

```cpp
// Get the diagonal block (local subdomain matrix)
mfem::HypreParMatrix* globalMatrix = ...;  // Your assembled parallel matrix
mfem::HypreParMatrix* localMatrix = globalMatrix->GetDiag();

// CRITICAL: localMatrix is a POINTER owned by globalMatrix
// DO NOT delete it manually!
```

**Memory Management Warning:**
- `GetDiag()` returns a pointer to internal data
- The returned matrix is owned by the parent `HypreParMatrix`
- Deleting it manually causes double-free errors
- It's automatically cleaned up when the parent is deleted

## One-Level Additive Schwarz Method

### Mathematical Formula

$$
M^{-1} = \sum_{i=1}^{N} R_i^T A_i^{-1} R_i
$$

Where:
- $N$ = number of subdomains (MPI ranks)
- $R_i$ = restriction operator to subdomain $i$
- $A_i$ = local subdomain matrix
- $A_i^{-1}$ = local subdomain solver (e.g., Gauss-Seidel)

### Implementation Pattern

```cpp
class OneLevelSchwarz : public mfem::Solver
{
private:
    mfem::HypreParMatrix* globalMatrix_;
    mfem::HypreSmoother* localSolver_;  // GS or Jacobi
    int myRank_;
    
public:
    OneLevelSchwarz(mfem::HypreParMatrix& A)
        : mfem::Solver(A.Height()), globalMatrix_(&A)
    {
        MPI_Comm_rank(A.GetComm(), &myRank_);
        
        // Extract local subdomain matrix
        mfem::HypreParMatrix* localA = globalMatrix_->GetDiag();
        
        // Create local solver (e.g., 3 sweeps of Gauss-Seidel)
        localSolver_ = new mfem::HypreSmoother(*localA, 
                                               mfem::HypreSmoother::GS, 3);
    }
    
    void Mult(const mfem::Vector& x, mfem::Vector& y) const override
    {
        // Apply local solve: y = A_i^{-1} * x
        localSolver_->Mult(x, y);
        
        // Communication is handled by the outer Krylov solver
        // via globalMatrix_->Mult()
    }
    
    ~OneLevelSchwarz()
    {
        delete localSolver_;
        // Do NOT delete globalMatrix_ (we don't own it)
    }
};
```

### Usage in CG Solver

```cpp
// Assemble global parallel matrix
mfem::HypreParMatrix A = ...;
mfem::Vector b = ...;
mfem::Vector x = ...;

// Create one-level Schwarz preconditioner
OneLevelSchwarz prec(A);

// Solve with PCG
mfem::CGSolver cg(MPI_COMM_WORLD);
cg.SetOperator(A);
cg.SetPreconditioner(prec);
cg.SetRelTol(1e-12);
cg.SetMaxIter(500);
cg.Mult(b, x);

std::cout << "CG iterations: " << cg.GetNumIterations() << std::endl;
```

### Performance Characteristics

**Problem:** Iteration count degrades with subdomain count.

| MPI Ranks | Subdomain Size | CG Iterations |
|-----------|----------------|---------------|
| 2         | Large          | 28            |
| 4         | Medium         | 32            |
| 8         | Small          | 47            |

**Why?** Information propagates only through subdomain overlaps. As subdomains shrink (more ranks), global modes take more iterations to resolve.

## Two-Level Schwarz Method (Framework)

### Mathematical Formula

$$
M^{-1} = \Psi A_H^{-1} \Psi^T + \sum_{i=1}^{N} R_i^T A_i^{-1} R_i
$$

**New term:** $\Psi A_H^{-1} \Psi^T$ = coarse-grid correction

- $\Psi$ = prolongation operator (coarse → fine DOFs)
- $A_H$ = coarse-grid matrix (via Galerkin: $A_H = \Psi^T A \Psi$)
- $A_H^{-1}$ = coarse-grid solver (AMG)

### Implementation Challenges in MFEM

**Challenge 1: Prolongation Operator Construction**

MFEM provides `ParFiniteElementSpace::Dof_TrueDof_Matrix()` for parallel DOF mappings, but building custom restriction/prolongation between different FE spaces is non-trivial:

```cpp
// Create coarse FE space (polynomial order 1)
mfem::H1_FECollection coarseFC(1, dim);
mfem::ParFiniteElementSpace coarseSpace(pmesh, &coarseFC);

// Get "built-in" prolongation from coarse to fine
mfem::HypreParMatrix* P = coarseSpace.Dof_TrueDof_Matrix();
// ^ This is NOT the prolongation we need!
// This is the local-to-true DOF mapping, not coarse-to-fine interpolation
```

**What's needed:** A custom interpolation operator that takes coefficients from coarse H1 space and evaluates them at fine H1 DOF locations.

**Challenge 2: Galerkin Projection**

Computing $A_H = \Psi^T A \Psi$ requires triple matrix product:

```cpp
// Pseudo-code (not actual MFEM API for custom operators)
mfem::HypreParMatrix* Psi_T = Psi->Transpose();
mfem::HypreParMatrix* temp = mfem::ParMult(A, Psi);
mfem::HypreParMatrix* A_H = mfem::ParMult(Psi_T, temp);
delete temp;
delete Psi_T;
```

This works when `Psi` is a proper `HypreParMatrix`, but constructing `Psi` itself is the bottleneck.

**Challenge 3: Memory Ownership**

```cpp
// INCORRECT - causes segfault:
~TwoLevelSchwarz()
{
    delete prolongation_;  // BAD! It's owned by coarseSpace_
    delete coarseSpace_;
}

// CORRECT:
~TwoLevelSchwarz()
{
    delete coarseSpace_;  // This cleans up prolongation_ too
}
```

### Current Framework Status

As of Phase 2.2, the two-level Schwarz class exists as a **framework**:

✅ **Implemented:**
- Class structure and interface
- Local subdomain solves (one-level part)
- Coarse FE space creation
- Coarse matrix assembly (simplified, not Galerkin)
- AMG solver for coarse problem

⚠️ **TODO (explicitly marked in code):**
- Proper prolongation operator $\Psi$ construction
- Galerkin projection $A_H = \Psi^T A \Psi$
- Restriction $\Psi^T x$ in `Mult()`
- Prolongation $\Psi y$ in `Mult()`

### Expected Performance (When Completed)

With proper coarse-grid correction, iteration count should remain constant:

| MPI Ranks | One-Level Iters | Two-Level Iters (Target) |
|-----------|-----------------|--------------------------|
| 2         | 28              | ~15-18                   |
| 4         | 32              | ~15-18                   |
| 8         | 47              | ~15-18                   |
| 16        | 80+             | ~15-18                   |

**Key property:** Coarse grid provides global communication channel, enabling $O(1)$ iteration count regardless of subdomain count.

## Parallel Communication Patterns

### Automatic Handling by MFEM

When you use `HypreParMatrix::Mult(x, y)`, MFEM automatically:

1. Packs shared DOFs into MPI buffers
2. Exchanges ghost values with neighboring ranks
3. Unpacks received data
4. Performs local SpMV
5. Reconciles boundary contributions

**You don't need to write MPI code!**

### Manual Communication (if needed)

```cpp
// Get parallel grid function
mfem::ParGridFunction u(&fespace);

// Exchange ghost values manually
u.ExchangeFaceNbrData();

// Now u has updated shared DOF values from neighbors
```

## Common Pitfalls

### Pitfall 1: Deleting Operator-Owned Pointers

```cpp
// BAD:
mfem::HypreParMatrix* diag = A.GetDiag();
delete diag;  // CRASH! A owns this pointer

// GOOD:
mfem::HypreParMatrix* diag = A.GetDiag();
// Use diag, but don't delete it
```

### Pitfall 2: Ignoring True vs. Local DOFs

```cpp
// BAD: Mixing local and true DOF vectors
mfem::Vector local_vec(fespace.GetVSize());  // Local size
A.Mult(local_vec, result);  // A expects true DOF size!

// GOOD: Use consistent true DOF sizes
mfem::Vector true_vec(fespace.GetTrueVSize());
A.Mult(true_vec, result);
```

### Pitfall 3: Incorrect Preconditioner Interface

```cpp
// Your preconditioner MUST inherit from mfem::Solver
class MyPrec : public mfem::Solver  // Correct base class
{
    void Mult(const mfem::Vector& x, mfem::Vector& y) const override
    {
        // Apply preconditioning: y = M^{-1} x
    }
    
    void SetOperator(const mfem::Operator& op) override
    {
        // Store reference to system matrix if needed
    }
};
```

## Performance Tips

### Tip 1: Tune Local Solver Iterations

```cpp
// More local smoothing sweeps = fewer global Krylov iterations
// But each sweep has cost, so balance is needed
localSolver_ = new mfem::HypreSmoother(*localA, 
                                       mfem::HypreSmoother::GS, 
                                       5);  // Try 1, 3, 5, 10
```

### Tip 2: Use AMG for Local Solves (if subdomains are large)

```cpp
// For large subdomains (>10k DOFs), AMG is better than GS
mfem::HypreBoomerAMG* localAMG = new mfem::HypreBoomerAMG(*localA);
localAMG->SetPrintLevel(0);
// Use localAMG instead of HypreSmoother
```

### Tip 3: Monitor Iteration Counts

```cpp
cg.SetPrintLevel(1);  // Print residual at each iteration

// After solve:
int iters = cg.GetNumIterations();
if (myRank == 0)
{
    std::cout << "Converged in " << iters << " iterations" << std::endl;
}
```

## References

- **MFEM Examples**: `ex26p.cpp` (parallel multigrid with HYPRE)
- **MFEM Miniapps**: `parelag` (element agglomeration for DDM)
- **Theory**: "Domain Decomposition Methods" (Smith, Bjørstad, Gropp, 1996)
- **Implementation**: `src/hpcfem/solver_one_level_schwarz.cpp`
- **Framework**: `src/hpcfem/solver_two_level_schwarz.cpp` (coarse correction TODO)

## Summary

**What Works:**
- One-level Schwarz with local Gauss-Seidel smoothers
- Automatic parallel mesh partitioning
- Local subdomain matrix extraction
- Integration with MFEM's parallel Krylov solvers

**What's Challenging:**
- Custom prolongation operator construction
- Galerkin projection for coarse operator
- Balancing local solve cost vs. global iterations

**Best Practice:**
For production code, use MFEM's built-in AMG (HYPRE BoomerAMG) or geometric multigrid unless you need specialized DDM features. DDM is most valuable for:
- Highly heterogeneous problems (domain-wise material properties)
- Problems where geometric coarsening is difficult
- Coupling with specialized subdomain solvers
