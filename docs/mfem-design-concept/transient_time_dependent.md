# Transient (Time-Dependent) Problems with MFEM

This guide covers solving time-dependent PDEs with MFEM, including parabolic (diffusion/heat) and hyperbolic (wave) equations.

## Overview

Transient problems involve PDEs with time derivatives. MFEM provides ODE time integration facilities to handle systems of the form:

$$\frac{du}{dt} = f(u,t)$$ (explicit form)

or

$$M \frac{du}{dt} = f(u,t)$$ (semi-discrete form)

where $M$ is the mass matrix.

## Problem Types

### 1. Parabolic (Diffusion/Heat Equation)

**Strong form:**
$$\frac{\partial u}{\partial t} - \nabla \cdot (k \nabla u) = f \quad \text{in } \Omega \times (0,T]$$

with initial condition $u(x,0) = u_0(x)$ and appropriate boundary conditions.

**Weak form (semi-discrete):**
$$M \frac{du}{dt} + K u = b$$

where:
- $M$ is the mass matrix: $(M)_{ij} = \int_\Omega \phi_i \phi_j$
- $K$ is the stiffness matrix: $(K)_{ij} = \int_\Omega k \nabla\phi_i \cdot \nabla\phi_j$
- $b$ is the load vector

### 2. Hyperbolic (Wave Equation)

**Strong form:**
$$\frac{\partial^2 u}{\partial t^2} - c^2 \nabla^2 u = f$$

Rewritten as first-order system:
$$\frac{\partial u}{\partial t} = v$$
$$\frac{\partial v}{\partial t} = c^2 \nabla^2 u + f$$

## MFEM Time Integration Framework

### Key Classes

1. **`TimeDependentOperator`**: Base class representing $\frac{du}{dt} = f(u,t)$
2. **`ODESolver`**: Abstract base for time integrators
3. **Concrete solvers**:
   - `ForwardEulerSolver`: Explicit Euler ($u^{n+1} = u^n + \Delta t f(u^n)$)
   - `RK4Solver`: 4th-order Runge-Kutta (explicit)
   - `BackwardEulerSolver`: Implicit Euler (requires solving nonlinear system)
   - `ImplicitMidpointSolver`: 2nd-order implicit
   - `SDIRK23Solver`, `SDIRK33Solver`: Implicit Runge-Kutta methods

### Time Integration Pattern

```cpp
// 1. Define the time-dependent operator
class HeatOperator : public TimeDependentOperator {
private:
    BilinearForm *M, *K;  // Mass and stiffness
    CGSolver M_solver;     // Inverse mass matrix solver
    
public:
    // Explicit: compute du/dt = M^{-1}(b - K*u)
    virtual void Mult(const Vector &u, Vector &dudt) const;
    
    // Implicit: solve (M + dt*K)*du = dt*rhs
    virtual void ImplicitSolve(const real_t dt, const Vector &u, Vector &du);
};

// 2. Setup initial condition
GridFunction u(&fes);
u.ProjectCoefficient(u0_coeff);

// 3. Choose time integrator
ODESolver *ode_solver;
ode_solver = new BackwardEulerSolver;  // or RK4Solver, etc.

// 4. Time loop
HeatOperator oper(M, K, ...);
ode_solver->Init(oper);
double t = 0.0, t_final = 1.0, dt = 0.01;
for (int ti = 0; t < t_final; ti++) {
    ode_solver->Step(u, t, dt);
    // Save/visualize u at time t
}
```

## Implementation Steps (Detailed)

### Step 1: Build Mass and Stiffness Forms

```cpp
// Mass matrix M
BilinearForm M(&fes);
ConstantCoefficient one(1.0);
M.AddDomainIntegrator(new MassIntegrator(one));
M.Assemble();
M.Finalize();

// Stiffness matrix K (diffusion)
BilinearForm K(&fes);
ConstantCoefficient k(1.0);  // conductivity
K.AddDomainIntegrator(new DiffusionIntegrator(k));
// Add boundary integrators if needed (Robin BC)
K.Assemble();
K.Finalize();
```

### Step 2: Define TimeDependentOperator

For heat equation: $M \frac{du}{dt} = -K u + b$

```cpp
class HeatOperator : public TimeDependentOperator {
private:
    SparseMatrix *M_mat, *K_mat;
    Vector b;  // RHS load
    CGSolver M_solver;
    DSmoother M_prec;
    mutable Vector z;
    
public:
    HeatOperator(BilinearForm &M_, BilinearForm &K_, const Vector &b_)
        : TimeDependentOperator(M_.Size()), b(b_)
    {
        M_mat = &M_.SpMat();
        K_mat = &K_.SpMat();
        z.SetSize(M_mat->Height());
        
        M_solver.SetOperator(*M_mat);
        M_solver.SetRelTol(1e-9);
        M_solver.SetMaxIter(100);
        M_solver.SetPrintLevel(0);
        M_solver.SetPreconditioner(M_prec);
    }
    
    // Explicit time stepping: du/dt = M^{-1}(-K*u + b)
    virtual void Mult(const Vector &u, Vector &dudt) const override {
        K_mat->Mult(u, z);     // z = K*u
        z.Neg();                // z = -K*u
        z += b;                 // z = -K*u + b
        M_solver.Mult(z, dudt); // dudt = M^{-1}*z
    }
    
    // Implicit solve: given dt, u_old, solve for k such that
    // u_new = u_old + dt*k, where M*k = -K*(u_old + dt*k) + b
    // Rearranging: (M + dt*K)*k = -K*u_old + b
    virtual void ImplicitSolve(const real_t dt, const Vector &u, Vector &k) override {
        // Form system matrix: A = M + dt*K
        SparseMatrix A = *M_mat;
        A.Add(dt, *K_mat);
        
        // RHS: rhs = -K*u + b
        Vector rhs(u.Size());
        K_mat->Mult(u, rhs);
        rhs.Neg();
        rhs += b;
        
        // Solve A*k = rhs
        CGSolver solver;
        DSmoother prec;
        solver.SetOperator(A);
        solver.SetRelTol(1e-9);
        solver.SetMaxIter(500);
        solver.SetPrintLevel(0);
        solver.SetPreconditioner(prec);
        solver.Mult(rhs, k);
    }
};
```

### Step 3: Time Loop with Visualization

```cpp
// Initial condition
FunctionCoefficient u0_func(u0_function);
GridFunction u(&fes);
u.ProjectCoefficient(u0_func);

// Create time-dependent operator
HeatOperator oper(M, K, b);

// Choose and initialize ODE solver
BackwardEulerSolver ode_solver;  // Or RK4Solver for explicit
ode_solver.Init(oper);

// Time parameters
double t = 0.0;
double dt = 0.01;
double t_final = 1.0;

// Output collection for time series
ParaViewDataCollection paraview("output_dir", &mesh);
paraview.SetLevelsOfDetail(order);
paraview.RegisterField("temperature", &u);
paraview.SetTime(t);
paraview.Save();

// Time loop
bool done = false;
for (int ti = 0; !done; ti++) {
    if (t + dt >= t_final - dt/2) {
        dt = t_final - t;
        done = true;
    }
    
    ode_solver.Step(u, t, dt);
    
    // Output at intervals
    if ((ti+1) % 10 == 0 || done) {
        std::cout << "time step " << ti+1 << ", t = " << t << std::endl;
        paraview.SetCycle(ti+1);
        paraview.SetTime(t);
        paraview.Save();
    }
}
```

## Boundary Conditions for Transient Problems

### Dirichlet (Essential)

Apply at each time step by projecting time-dependent coefficient:

```cpp
// Time-dependent BC: u = g(x,t) on Gamma_D
class TimeDependentBC : public Coefficient {
    double t_current;
public:
    void SetTime(double t) { t_current = t; }
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
        Vector x; T.Transform(ip, x);
        return g(x[0], x[1], t_current);  // User-defined function
    }
};

// In time loop:
TimeDependentBC dirichlet_bc;
for (int ti = 0; !done; ti++) {
    dirichlet_bc.SetTime(t);
    u.ProjectBdrCoefficient(dirichlet_bc, ess_bdr);
    ode_solver.Step(u, t, dt);
}
```

### Neumann (Natural)

Include in RHS vector $b$ which may be time-dependent:

```cpp
// Update load vector at each time step
LinearForm b(&fes);
TimeDependentNeumann flux_bc;
for (int ti = 0; !done; ti++) {
    flux_bc.SetTime(t);
    b = 0.0;
    b.AddBoundaryIntegrator(new BoundaryLFIntegrator(flux_bc), neumann_bdr);
    b.Assemble();
    // Use updated b in oper
    ode_solver.Step(u, t, dt);
}
```

### Robin (Mixed)

Add to stiffness matrix (time-independent coefficient) or update at each step.

## Stability and Time Step Selection

### Explicit Methods (Forward Euler, RK)

**Stability constraint (CFL condition):**
$$\Delta t \leq C \frac{h^2}{k}$$

where $h$ is mesh size, $k$ is diffusivity, and $C$ is a constant (typically $C \approx 0.5$ for Forward Euler).

**Use cases:**
- Small problems
- Problems dominated by explicit source terms
- Hyperbolic equations with appropriate CFL

### Implicit Methods (Backward Euler, SDIRK)

**Unconditionally stable** for linear diffusion (no CFL constraint).

**Accuracy constraint:**
$$\Delta t \sim O(h^p)$$
where $p$ is the spatial discretization order for accuracy.

**Use cases:**
- Stiff problems (large diffusivity)
- Long-time integration
- Large time steps needed

## Performance Optimization

### 1. Mass Matrix Lumping (Explicit only)

Replace $M$ with diagonal approximation $M_L$:

```cpp
BilinearForm M(&fes);
M.AddDomainIntegrator(new MassIntegrator);
M.Assemble();
M.Finalize();

// Lump (replace with diagonal)
Vector diag;
M.SpMat().GetDiag(diag);
M_lumped = new SparseMatrix;  // Diagonal matrix
```

Lumped mass matrix allows explicit inversion without solving linear system.

### 2. Partial Assembly

For matrix-free operations:

```cpp
BilinearForm K(&fes);
K.SetAssemblyLevel(AssemblyLevel::PARTIAL);
K.AddDomainIntegrator(new DiffusionIntegrator(k));
K.Assemble();
// K can be applied without forming sparse matrix
```

### 3. Adaptive Time Stepping

Use error estimators to adapt $\Delta t$:

```cpp
// Use SDIRK or RK methods with error estimation
SDIRK23Solver ode_solver;  // 2nd order with 3rd order error estimate
ode_solver.Init(oper);

// Set tolerances
double reltol = 1e-4, abstol = 1e-6;
// Adaptive stepping handled internally
```

## Common Pitfalls

1. **Forgetting to assemble and finalize forms** before creating operators
2. **Not updating time-dependent BCs** at each time step
3. **Stability issues with explicit methods**: reduce time step
4. **Mass matrix singularity**: ensure proper BCs are applied
5. **Performance**: for large 3D problems, use implicit methods with good preconditioners (AMG)

## Example Mapping

- **Explicit heat equation**: See `example/ex3_transient_heat.cpp` (to be created)
- **MFEM official examples**:
  - `ex9.cpp`: DG advection (explicit time stepping)
  - `ex10.cpp`: Nonlinear elastodynamics (implicit)
  - `ex16.cpp`: Time-dependent Navier-Stokes

## Summary Checklist

- [ ] Build mass matrix $M$ and stiffness matrix $K$
- [ ] Define initial condition $u_0$
- [ ] Create `TimeDependentOperator` subclass implementing `Mult` and/or `ImplicitSolve`
- [ ] Choose appropriate `ODESolver` (explicit vs implicit based on stability)
- [ ] Set up time loop with proper time step
- [ ] Handle time-dependent boundary conditions
- [ ] Use `ParaViewDataCollection` for time series output
- [ ] Verify stability and accuracy by varying $\Delta t$ and mesh resolution

## References

- MFEM examples: `ex9.cpp`, `ex10.cpp`, `ex16.cpp`
- MFEM documentation: Time-dependent operators and ODE integrators
- Theory: Finite Element Methods for PDEs, Quarteroni & Valli
