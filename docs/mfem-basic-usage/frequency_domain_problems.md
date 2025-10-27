# Frequency-Domain Problems with MFEM

This guide covers solving frequency-domain (harmonic) PDE problems with MFEM, particularly electromagnetics and acoustics applications.

## Overview

Frequency-domain problems arise when time-harmonic solutions $u(x,t) = \tilde{u}(x) e^{j\omega t}$ are sought, where $\omega = 2\pi f$ is angular frequency. The time derivative becomes a multiplication: $\frac{\partial}{\partial t} \to j\omega$.

## Problem Types

### 1. Helmholtz Equation (Acoustics, EM Scalar)

**Strong form:**
$$\nabla^2 u + k^2 u = f \quad \text{in } \Omega$$

where $k = \omega/c$ is the wavenumber ($c$ = wave speed).

**Weak form:**
$$\int_\Omega (\nabla u \cdot \nabla v - k^2 u v) \, dx = \int_\Omega f v \, dx + \text{BC terms}$$

**Discrete form:**
$$(-K + k^2 M) u = b$$

**Applications:**

- Acoustic wave propagation
- Scalar electromagnetic fields
- Quantum scattering

### 2. Vector Helmholtz / Time-Harmonic Maxwell

**Strong form (electric field formulation):**
$$\nabla \times (\mu_r^{-1} \nabla \times E) - k_0^2 \epsilon_r E = -j\omega J$$

where $k_0 = \omega/c_0$ is free-space wavenumber.

**Weak form:**
$$\int_\Omega \mu_r^{-1} (\nabla \times E) \cdot (\nabla \times F) \, dx - k_0^2 \int_\Omega \epsilon_r E \cdot F \, dx = \int_\Omega (-j\omega J) \cdot F \, dx + \text{BC}$$

**Discrete form:**
$$(C - k_0^2 M_\epsilon) E = b$$

where $C$ is curl-curl matrix, $M_\epsilon$ is mass matrix with permittivity.

**Applications:**

- Waveguides and cavities
- Antennas and scattering
- Photonics and metamaterials

### 3. Driven Resonator / Port Excitation

Modal excitation at boundaries (ports) with impedance matching or absorbing BCs.

## MFEM Implementation Patterns

### Helmholtz Equation (Scalar)

```cpp
// 1. Setup
Mesh mesh("mesh.mesh");
int order = 2;
H1_FECollection fec(order, mesh.Dimension());
FiniteElementSpace fes(&mesh, &fec);

// 2. Frequency parameters
double freq = 1e9;  // 1 GHz
double c0 = 299792458.0;  // speed of light [m/s]
double k0 = 2.0 * M_PI * freq / c0;  // wavenumber

// 3. Bilinear form: -grad·grad + k^2 * mass
BilinearForm a(&fes);
ConstantCoefficient one(1.0);
a.AddDomainIntegrator(new DiffusionIntegrator(one));  // -nabla^2 term
ConstantCoefficient k_sq(k0 * k0);
a.AddDomainIntegrator(new MassIntegrator(k_sq));  // k^2 term
a.Assemble();
a.Finalize();

// 4. Linear form (RHS)
LinearForm b(&fes);
// Source term or boundary excitation
FunctionCoefficient source_func(source_function);
b.AddDomainIntegrator(new DomainLFIntegrator(source_func));
b.Assemble();

// 5. Boundary conditions
Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 0;
ess_bdr[port_attr-1] = 1;  // Driven port
Array<int> ess_tdof_list;
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

GridFunction u(&fes);
u = 0.0;
FunctionCoefficient port_excitation(port_function);
u.ProjectBdrCoefficient(port_excitation, ess_bdr);

// 6. Form linear system
SparseMatrix A;
Vector B, X;
a.FormLinearSystem(ess_tdof_list, u, b, A, X, B);

// 7. Solve (indefinite system, use direct or GMRES)
#ifdef MFEM_USE_SUITESPARSE
UMFPackSolver umf;
umf.SetOperator(A);
umf.Mult(B, X);
#else
GMRESSolver gmres;
GSSmoother prec(A);
gmres.SetOperator(A);
gmres.SetRelTol(1e-8);
gmres.SetMaxIter(500);
gmres.SetPrintLevel(1);
gmres.SetPreconditioner(prec);
gmres.Mult(B, X);
#endif

// 8. Recover solution
a.RecoverFEMSolution(X, b, u);
```

### Vector Helmholtz (Electromagnetics)

```cpp
// 1. Setup with Nedelec (H(curl)) elements
Mesh mesh("waveguide.mesh");
int order = 1;
ND_FECollection fec(order, mesh.Dimension());
FiniteElementSpace fes(&mesh, &fec);

// 2. Parameters
double freq = 10e9;  // 10 GHz
double c0 = 299792458.0;
double omega = 2.0 * M_PI * freq;
double k0 = omega / c0;
double eps_r = 1.0;  // Relative permittivity
double mu_r = 1.0;   // Relative permeability

// 3. Bilinear form: curl-curl - k0^2 * eps_r * mass
BilinearForm a(&fes);
ConstantCoefficient one_over_mu(1.0 / mu_r);
a.AddDomainIntegrator(new CurlCurlIntegrator(one_over_mu));

double mass_coef = -(k0 * k0 * eps_r);  // Note: negative sign
ConstantCoefficient neg_k_sq_eps(mass_coef);
a.AddDomainIntegrator(new VectorFEMassIntegrator(neg_k_sq_eps));
a.Assemble();
a.Finalize();

// 4. Essential BCs (PEC walls and ports)
Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 0;
// PEC (Perfect Electric Conductor): n x E = 0
ess_bdr[pec_attr_1-1] = 1;
ess_bdr[pec_attr_2-1] = 1;
// Input port
ess_bdr[port_in-1] = 1;
// Output port (could be zero or matched load)
ess_bdr[port_out-1] = 1;

Array<int> ess_tdof_list;
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

// 5. Project port excitation (modal field)
GridFunction E(&fes);
E = 0.0;

// Input port: project modal field (e.g., TE10 mode)
Array<int> port_in_bdr(mesh.bdr_attributes.Max());
port_in_bdr = 0;
port_in_bdr[port_in-1] = 1;
VectorFunctionCoefficient te10_mode(3, te10_modal_field);
E.ProjectBdrCoefficientTangent(te10_mode, port_in_bdr);

// Output port: zero or matched (absorbing)
Array<int> port_out_bdr(mesh.bdr_attributes.Max());
port_out_bdr = 0;
port_out_bdr[port_out-1] = 1;
Vector zero_vec(3); zero_vec = 0.0;
VectorConstantCoefficient zero_field(zero_vec);
E.ProjectBdrCoefficientTangent(zero_field, port_out_bdr);

// 6. Form system and solve
LinearForm rhs(&fes);
rhs = 0.0;
rhs.Assemble();

SparseMatrix A;
Vector B, X;
a.FormLinearSystem(ess_tdof_list, E, rhs, A, X, B);

// Solve with direct solver (preferred for indefinite systems)
#ifdef MFEM_USE_SUITESPARSE
UMFPackSolver umf;
umf.SetOperator(A);
umf.Mult(B, X);
#else
// Fallback to GMRES
GMRESSolver gmres;
GSSmoother prec(A);
gmres.SetOperator(A);
gmres.SetRelTol(1e-8);
gmres.SetMaxIter(1000);
gmres.SetPrintLevel(1);
gmres.SetPreconditioner(prec);
gmres.Mult(B, X);
#endif

a.RecoverFEMSolution(X, rhs, E);
```

## Boundary Conditions

### 1. Dirichlet (Essential)

**Scalar:** $u = g$ on $\Gamma_D$
**Vector:** $n \times E = g$ on $\Gamma_D$ (tangential component prescribed)

Implementation: project using `ProjectBdrCoefficient` (scalar) or `ProjectBdrCoefficientTangent` (vector).

### 2. Neumann (Natural)

**Scalar:** $\nabla u \cdot n = h$ on $\Gamma_N$
**Vector:** $n \times (\nabla \times E) = h$ on $\Gamma_N$

Implementation: add to RHS via `BoundaryLFIntegrator`.

### 3. Robin / Impedance BC

**Form:** $\nabla u \cdot n + \alpha u = g$ on $\Gamma_R$

Weak form contribution:

$$\int_{\Gamma_R} \alpha u v \, ds = \int_{\Gamma_R} g v \, ds$$

Implementation:

```cpp
// Add to bilinear form (LHS)
ConstantCoefficient alpha(impedance_value);
Array<int> robin_bdr(mesh.bdr_attributes.Max());
robin_bdr = 0; robin_bdr[robin_attr-1] = 1;
a.AddBoundaryIntegrator(new MassIntegrator(alpha), robin_bdr);

// Add to linear form (RHS)
ConstantCoefficient g(excitation_value);
b.AddBoundaryIntegrator(new BoundaryLFIntegrator(g), robin_bdr);
```

### 4. Absorbing BC (Sommerfeld / PML)

For open-domain problems to prevent reflections.

**First-order absorbing:**
$$\nabla u \cdot n + jk u = 0$$

Implemented as Robin BC with complex coefficient $\alpha = jk$.

**Perfectly Matched Layer (PML):**
Use complex coordinate stretching in boundary layer. MFEM supports complex arithmetic:

```cpp
// Enable complex arithmetic
#define MFEM_USE_COMPLEX_OPERATOR
// Define complex-valued problem
```

## Complex-Valued Problems

Many frequency-domain problems are inherently complex-valued. MFEM handles this in two ways:

### Method 1: Real-Imaginary Splitting

Split $u = u_r + j u_i$ into real and imaginary parts, forming a $2N \times 2N$ real system:

$$\begin{bmatrix} A & -B \\ B & A \end{bmatrix} \begin{bmatrix} u_r \\ u_i \end{bmatrix} = \begin{bmatrix} f_r \\ f_i \end{bmatrix}$$

(Not commonly used in MFEM; complex problems often solved with native complex support externally or via block systems.)

### Method 2: External Complex Solvers

Use PETSc or other libraries with native complex support:

```cpp
#ifdef MFEM_USE_PETSC
// PETSc can handle complex-valued problems natively
PetscLinearSolver petsc_solver(MPI_COMM_WORLD);
// ... configure for complex system
#endif
```

### Method 3: Magnitude-Phase Formulation

For lossless problems, sometimes sufficient to solve real-valued system at resonance and compute phase separately. (Application-specific.)

## Modal Analysis of Waveguides

For waveguide problems, often interested in propagation constant $\beta$ (eigenvalue) for given frequency:

$$(\nabla_t^2 + k^2 - \beta^2) E_t = 0$$

This becomes a 2D eigenvalue problem on the cross-section.

**Implementation:**

- Solve 2D eigenvalue problem: $(K_{t} - \beta^2 M) E_t = 0$
- For each mode $i$, $\beta_i$ is the propagation constant
- Full 3D field: $E(x,y,z) = E_t(x,y) e^{-j\beta z}$

## Performance Considerations

### 1. Direct vs Iterative Solvers

**Direct solvers (UMFPACK, MUMPS, SuperLU):**

- Reliable for indefinite systems
- Good for moderate problem sizes ($N < 10^6$)
- High memory usage: $O(N^{1.5})$ to $O(N^2)$

**Iterative solvers (GMRES, BiCGSTAB):**

- Lower memory: $O(N)$
- Require good preconditioner
- Indefinite Helmholtz is challenging (large $k$)

### 2. Preconditioners for Helmholtz

**Challenges:**

- Standard AMG fails for indefinite problems
- Complex shifted Laplacian: $(K + j\alpha k^2 M)^{-1}$ as preconditioner

**Strategies:**

```cpp
// Shifted Laplacian preconditioner
BilinearForm prec_form(&fes);
ConstantCoefficient one(1.0);
prec_form.AddDomainIntegrator(new DiffusionIntegrator(one));
double shift = 0.5;  // Tuning parameter
ConstantCoefficient shift_mass(shift * k0 * k0);
prec_form.AddDomainIntegrator(new MassIntegrator(shift_mass));
prec_form.Assemble();
prec_form.Finalize();

// Use AMG on shifted operator
HypreBoomerAMG prec(prec_form.SpMat());
gmres.SetPreconditioner(prec);
```

### 3. Matrix-Free / Partial Assembly

For large 3D problems:

```cpp
BilinearForm a(&fes);
a.SetAssemblyLevel(AssemblyLevel::PARTIAL);
a.AddDomainIntegrator(new CurlCurlIntegrator(one_over_mu));
a.AddDomainIntegrator(new VectorFEMassIntegrator(neg_k_sq_eps));
a.Assemble();
// Matrix-free matvec via operator application
```

## Validation and Post-Processing

### S-Parameters (Scattering Parameters)

For waveguide/cavity problems:

$$S_{21} = \frac{\text{Power out (port 2)}}{\text{Power in (port 1)}}$$

Compute from field integrals at ports.

### Quality Factor (Q)

For resonant cavities:
$$Q = \frac{\text{Stored energy}}{\text{Power loss}}$$

Computed from field distribution.

### Field Visualization

```cpp
ParaViewDataCollection paraview("solution", &mesh);
paraview.RegisterField("E_field", &E);
paraview.Save();
```

For complex fields, save magnitude and phase separately:

```cpp
// Magnitude: |E|
GridFunction E_mag(&fes);
for (int i = 0; i < E.Size(); i++) {
    E_mag[i] = std::abs(E[i]);  // If E is complex-valued
}
paraview.RegisterField("E_magnitude", &E_mag);
```

## Common Pitfalls

1. **Sign convention:** Different authors use $e^{+j\omega t}$ vs $e^{-j\omega t}$. Affects sign of imaginary parts.
2. **Indefinite systems:** Standard preconditioners (AMG) don't work well. Use direct solver or specialized preconditioners.
3. **Spurious modes:** In H(curl) problems, ensure discrete kernel is handled (penalty on $\nabla \cdot E$ or tree-cotree gauge).
4. **High frequency:** Large $k$ makes problem increasingly indefinite and harder to solve iteratively.
5. **Boundary reflections:** Without proper absorbing BC or PML, artificial reflections contaminate solution.

## Example Mapping

- **Rectangular waveguide**: See `example/ex2_rect_waveguide.cpp` (existing)
- **MFEM official examples:**
  - `ex13p.cpp`: Maxwell cavity eigenvalues (related to frequency domain)
  - `ex22.cpp`: Complex-valued Helmholtz

## Summary Checklist

- [ ] Define frequency/wavenumber parameters ($k_0 = \omega/c_0$)
- [ ] Choose appropriate finite element space (H1 for scalar, H(curl) for vector)
- [ ] Build curl-curl or Laplacian + mass term with correct signs
- [ ] Apply boundary conditions (PEC, ports, absorbing)
- [ ] Project modal excitation at driven ports
- [ ] Solve indefinite system (prefer direct solver for reliability)
- [ ] Post-process: compute S-parameters, Q-factor, field distributions
- [ ] Visualize magnitude and phase of solution

## References

- MFEM examples: `ex22.cpp`
- Jin, J.M., "The Finite Element Method in Electromagnetics"
- Monk, P., "Finite Element Methods for Maxwell's Equations"
- Ihlenburg, F., "Finite Element Analysis of Acoustic Scattering"
