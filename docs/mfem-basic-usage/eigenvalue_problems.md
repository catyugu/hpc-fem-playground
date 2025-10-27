# Eigenvalue Problems with MFEM

This guide covers solving eigenvalue problems (EVPs) with MFEM, focusing on generalized eigenvalue problems arising from PDE discretizations.

## Overview

Eigenvalue problems in finite element analysis have the form:

$$K u = \lambda M u$$

where:

- $K$ is the stiffness matrix (often represents $-\nabla^2$ or similar operator)
- $M$ is the mass matrix (identity-like, represents $L^2$ inner product)
- $\lambda$ are eigenvalues
- $u$ are eigenvectors (eigenmodes/eigenfunctions)

## Problem Types

### 1. Laplacian Eigenvalue Problem

**Strong form:**
$$-\nabla^2 u = \lambda u \quad \text{in } \Omega$$

with homogeneous Dirichlet boundary conditions: $u = 0$ on $\partial\Omega$.

**Applications:**

- Vibration modes of a membrane
- Quantum mechanics (particle in a box)
- Heat equation characteristic functions

**Weak form:**
$$\int_\Omega \nabla u \cdot \nabla v \, dx = \lambda \int_\Omega u v \, dx$$

Discrete form:
$$K u = \lambda M u$$

### 2. Maxwell Eigenvalue Problem (Electromagnetic Cavities)

**Strong form:**
$$\nabla \times (\nabla \times E) = \lambda \epsilon E$$

with boundary condition $n \times E = 0$ (PEC walls).

**Applications:**

- Resonant modes in RF cavities
- Photonic crystals

**Discrete form (H(curl) space):**
$$C u = \lambda M u$$

where $C$ is curl-curl matrix, $M$ is mass matrix in $H(\text{curl})$.

### 3. Elastic Vibration Modes

**Strong form:**
$$-\nabla \cdot \sigma(u) = \lambda \rho u$$

where $\sigma$ is stress tensor, $\rho$ is density.

**Applications:**

- Structural vibration analysis
- Modal analysis in mechanical engineering

## MFEM Eigenvalue Solver Framework

### Key Components

1. **Eigenvalue Solvers:**
   - `HypreLOBPCG`: Locally Optimal Block Preconditioned Conjugate Gradient (requires HYPRE)
   - `HypreAME`: Auxiliary space Maxwell Eigensolver (for H(curl) problems)
   - `SLEPcEigenSolver`: Interface to SLEPc (requires PETSc/SLEPc)

2. **Preconditioners:**
   - `HypreBoomerAMG`: Algebraic multigrid for H1 problems
   - `HypreAMS`: Auxiliary space for H(curl) problems

### Standard Workflow

```cpp
// 1. Build stiffness matrix K
BilinearForm K(&fes);
K.AddDomainIntegrator(new DiffusionIntegrator);
K.Assemble();
K.Finalize();

// 2. Build mass matrix M
BilinearForm M(&fes);
M.AddDomainIntegrator(new MassIntegrator);
M.Assemble();
M.Finalize();

// 3. Handle essential BCs (eliminate rows/cols)
Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 1;  // Homogeneous Dirichlet on all boundaries
Array<int> ess_tdof_list;
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

// Form constrained operators
HypreParMatrix *K_mat, *M_mat;
K.FormSystemMatrix(ess_tdof_list, K_mat);
M.FormSystemMatrix(ess_tdof_list, M_mat);

// 4. Setup eigensolver
HypreBoomerAMG *amg = new HypreBoomerAMG(*K_mat);
amg->SetPrintLevel(0);

HypreLOBPCG *lobpcg = new HypreLOBPCG(MPI_COMM_WORLD);
lobpcg->SetNumModes(nev);  // Number of eigenvalues to compute
lobpcg->SetPreconditioner(*amg);
lobpcg->SetMaxIter(200);
lobpcg->SetTol(1e-8);
lobpcg->SetPrecondUsageMode(1);
lobpcg->SetPrintLevel(1);
lobpcg->SetMassMatrix(*M_mat);
lobpcg->SetOperator(*K_mat);

// 5. Set random initial vectors
lobpcg->SetInitialVectors(nev, seed);

// 6. Solve
lobpcg->Solve();

// 7. Extract eigenvalues and eigenvectors
Array<double> eigenvalues;
lobpcg->GetEigenvalues(eigenvalues);

for (int i = 0; i < nev; i++) {
    GridFunction u(&fes);
    u = lobpcg->GetEigenvector(i);
    // Visualize or save u
}
```

## Implementation Details

### Step 1: Finite Element Space

For Laplacian eigenvalue problem, use H1 space:

```cpp
Mesh mesh("mesh.mesh");
int order = 1;  // Polynomial order
H1_FECollection fec(order, mesh.Dimension());
FiniteElementSpace fes(&mesh, &fec);
```

For Maxwell, use Nedelec (H(curl)) space:

```cpp
ND_FECollection fec(order, mesh.Dimension());
FiniteElementSpace fes(&mesh, &fec);
```

### Step 2: Bilinear Forms

**Laplacian problem:**

```cpp
BilinearForm K(&fes);
ConstantCoefficient one(1.0);
K.AddDomainIntegrator(new DiffusionIntegrator(one));
K.Assemble();
K.Finalize();

BilinearForm M(&fes);
M.AddDomainIntegrator(new MassIntegrator(one));
M.Assemble();
M.Finalize();
```

**Maxwell problem:**

```cpp
BilinearForm C(&fes);
ConstantCoefficient one(1.0);
C.AddDomainIntegrator(new CurlCurlIntegrator(one));
C.Assemble();
C.Finalize();

BilinearForm M(&fes);
ConstantCoefficient epsilon(1.0);  // Permittivity
M.AddDomainIntegrator(new VectorFEMassIntegrator(epsilon));
M.Assemble();
M.Finalize();
```

### Step 3: Boundary Conditions

For homogeneous Dirichlet (most common in EVPs):

```cpp
Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 1;  // All boundaries
Array<int> ess_tdof_list;
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

// Eliminate constrained DOFs from matrices
HypreParMatrix *K_mat, *M_mat;
K.FormSystemMatrix(ess_tdof_list, K_mat);
M.FormSystemMatrix(ess_tdof_list, M_mat);
```

For periodic BCs, use periodic mesh (no explicit BC marking).

### Step 4: Solving with LOBPCG

```cpp
// Preconditioner (critical for performance)
HypreBoomerAMG *amg = new HypreBoomerAMG(*K_mat);
amg->SetPrintLevel(0);

// LOBPCG solver
HypreLOBPCG *lobpcg = new HypreLOBPCG(MPI_COMM_WORLD);
lobpcg->SetNumModes(nev);  // e.g., 10 smallest eigenvalues
lobpcg->SetPreconditioner(*amg);
lobpcg->SetMaxIter(200);
lobpcg->SetTol(1e-8);
lobpcg->SetPrecondUsageMode(1);
lobpcg->SetPrintLevel(1);

// Set operators
lobpcg->SetMassMatrix(*M_mat);
lobpcg->SetOperator(*K_mat);

// Random seed for initial guess
lobpcg->SetInitialVectors(nev, seed);

// Solve K u = lambda M u
lobpcg->Solve();

// Extract results
Array<double> eigenvalues;
lobpcg->GetEigenvalues(eigenvalues);

for (int i = 0; i < nev; i++) {
    std::cout << "Eigenvalue " << i << ": " << eigenvalues[i] << std::endl;
}
```

### Step 5: Extract and Visualize Eigenmodes

```cpp
ParaViewDataCollection paraview("eigenmodes", &mesh);
paraview.SetLevelsOfDetail(order);

for (int i = 0; i < nev; i++) {
    GridFunction u(&fes);
    u = lobpcg->GetEigenvector(i);
    
    // Normalize
    double norm = u.Norml2();
    u /= norm;
    
    // Save
    std::ostringstream mode_name;
    mode_name << "mode_" << i;
    paraview.RegisterField(mode_name.str().c_str(), &u);
    paraview.SetCycle(i);
    paraview.SetTime(eigenvalues[i]);  // Store eigenvalue as "time"
    paraview.Save();
}
```

## Theoretical Background

### Rayleigh Quotient

For a given vector $v$, the Rayleigh quotient is:

$$R(v) = \frac{v^T K v}{v^T M v}$$

**Properties:**

- $R(u_i) = \lambda_i$ for eigenvector $u_i$
- $\lambda_{\min} \leq R(v) \leq \lambda_{\max}$ for any $v$
- Minimizing $R(v)$ gives smallest eigenvalue

### LOBPCG Algorithm

**Locally Optimal Block Preconditioned Conjugate Gradient:**

1. Start with random initial block of vectors
2. Iterate:
   - Apply preconditioner
   - Compute Rayleigh-Ritz approximation in subspace
   - Update search directions (conjugate gradient style)
3. Converge to smallest eigenvalues

**Advantages:**

- Finds multiple eigenvalues simultaneously
- Good preconditioner → fast convergence
- Memory efficient (no full matrix factorization)

### Convergence

Number of iterations depends on:

- Spectral gap: $(\lambda_{k+1} - \lambda_k) / \lambda_k$
- Preconditioner quality
- Initial guess quality

Typical convergence: $10-100$ iterations for small eigenvalues with good AMG preconditioner.

## Performance Optimization

### 1. Use Good Preconditioner

For H1 (Laplacian):

```cpp
HypreBoomerAMG *amg = new HypreBoomerAMG(*K_mat);
amg->SetPrintLevel(0);
// Tuning for eigenvalue problems:
amg->SetCoarsening(10);  // HMIS coarsening
amg->SetInterpolation(6);  // Extended+i
amg->SetRelaxType(8);  // l1-Gauss-Seidel
```

For H(curl) (Maxwell):

```cpp
HypreAMS *ams = new HypreAMS(*C_mat, &fes);
ams->SetPrintLevel(0);
```

### 2. Partial Assembly (Matrix-Free)

For very large problems:

```cpp
BilinearForm K(&fes);
K.SetAssemblyLevel(AssemblyLevel::PARTIAL);
K.AddDomainIntegrator(new DiffusionIntegrator);
K.Assemble();
// Use matrix-free operators in eigensolver
```

### 3. Parallel Execution

Use parallel MFEM and MPI:

```cpp
ParMesh pmesh(MPI_COMM_WORLD, mesh);
ParFiniteElementSpace fes(&pmesh, &fec);
// All subsequent operations are parallel
```

## Common Pitfalls

1. **Spurious zero eigenvalues**: Check that essential BCs are properly applied
2. **Poor convergence**: Improve preconditioner settings
3. **Wrong eigenvalues**: Ensure problem is well-posed (check mesh quality, BCs)
4. **Maxwell spurious modes**: For H(curl), zero-frequency modes may appear; filter by checking $\|\nabla \cdot E\|$

## Special Cases

### Buckling Analysis

Problem: $K u = \lambda K_G u$ where $K_G$ is geometric stiffness.

Replace $M$ with $K_G$ in standard workflow.

### Shift-Invert Spectral Transformation

To find eigenvalues near a target $\sigma$:

Solve $(K - \sigma M)^{-1} M u = \mu u$

where eigenvalues are $\lambda = \sigma + 1/\mu$.

(Requires factorization or iterative solver for $(K - \sigma M)^{-1}$.)

## Validation

For simple geometries (square, cube, sphere), analytical solutions exist:

**2D square $[0, \pi] \times [0, \pi]$ with Dirichlet BC:**
$$\lambda_{m,n} = m^2 + n^2, \quad m,n = 1,2,3,...$$

**3D cube $[0, \pi]^3$ with Dirichlet BC:**
$$\lambda_{l,m,n} = l^2 + m^2 + n^2, \quad l,m,n = 1,2,3,...$$

Check computed eigenvalues against these for verification.

## Example Mapping

- **Laplacian eigenvalue**: See `example/ex4_eigenvalue_laplacian.cpp` (to be created)
- **MFEM official examples**:
  - `ex11p.cpp`: Parallel Laplacian eigenvalue problem
  - `ex13p.cpp`: Maxwell eigenvalue problem (cavity modes)

## Summary Checklist

- [ ] Build stiffness matrix $K$ and mass matrix $M$
- [ ] Apply essential boundary conditions
- [ ] Choose appropriate preconditioner (AMG for H1, AMS for H(curl))
- [ ] Configure LOBPCG solver (number of modes, tolerance)
- [ ] Set random initial vectors
- [ ] Solve and extract eigenvalues/eigenvectors
- [ ] Validate against analytical solutions (if available)
- [ ] Visualize eigenmodes in ParaView

## References

- MFEM examples: `ex11p.cpp`, `ex13p.cpp`
- HYPRE documentation: LOBPCG solver
- Knyazev, A.V., "Toward the Optimal Preconditioned Eigensolver: Locally Optimal Block Preconditioned Conjugate Gradient Method"
