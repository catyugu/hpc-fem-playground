# Comprehensive Boundary Conditions Guide for MFEM

This document provides a systematic guide to implementing various boundary conditions in MFEM across different problem types and finite element spaces.

## Overview

MFEM classifies boundary conditions into two main categories:

1. **Essential (Dirichlet)**: Prescribed solution values - enforced strongly by constraining DOFs
2. **Natural (Neumann, Robin)**: Prescribed fluxes or mixed conditions - enforced weakly through boundary integrals

## Boundary Attribute Marking

MFEM uses integer **boundary attributes** to mark different boundary regions in the mesh.

### Checking Available Attributes

```cpp
Mesh mesh("mesh.mesh");
std::cout << "Number of boundary attributes: " << mesh.bdr_attributes.Max() << std::endl;
for (int i = 0; i < mesh.bdr_attributes.Size(); i++) {
    std::cout << "Attribute: " << mesh.bdr_attributes[i] << std::endl;
}
```

### Creating Attribute Arrays

```cpp
// Array indexed 0 to (max_attr - 1), but refers to attributes 1 to max_attr
Array<int> bdr_marker(mesh.bdr_attributes.Max());
bdr_marker = 0;  // Initialize all to "not marked"
bdr_marker[attr_number - 1] = 1;  // Mark attribute attr_number
```

**Important:** Attributes are 1-indexed in mesh files, but arrays are 0-indexed in MFEM.

## Essential (Dirichlet) Boundary Conditions

### Scalar Problems (H1 Space)

#### Homogeneous: $u = 0$ on $\Gamma_D$

```cpp
Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 0;
ess_bdr[attr - 1] = 1;  // Mark boundary

Array<int> ess_tdof_list;
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

GridFunction u(&fes);
u = 0.0;  // Homogeneous value

// Form system with constraints
SparseMatrix A;
Vector B, X;
a.FormLinearSystem(ess_tdof_list, u, b, A, X, B);
```

#### Non-homogeneous: $u = g(x)$ on $\Gamma_D$

```cpp
// Define boundary value function
double dirichlet_function(const Vector &x) {
    return 300.0;  // Temperature, for example
}

FunctionCoefficient g(dirichlet_function);

Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 0;
ess_bdr[attr - 1] = 1;

Array<int> ess_tdof_list;
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

GridFunction u(&fes);
u = 0.0;
u.ProjectBdrCoefficient(g, ess_bdr);  // Project g onto boundary DOFs

// Form system (u already contains BC values)
a.FormLinearSystem(ess_tdof_list, u, b, A, X, B);
```

#### Multiple Regions with Different Values

```cpp
ConstantCoefficient cold(273.15);
ConstantCoefficient hot(373.15);

Array<int> cold_bdr(mesh.bdr_attributes.Max());
cold_bdr = 0; cold_bdr[1-1] = 1;
u.ProjectBdrCoefficient(cold, cold_bdr);

Array<int> hot_bdr(mesh.bdr_attributes.Max());
hot_bdr = 0; hot_bdr[3-1] = 1;
u.ProjectBdrCoefficient(hot, hot_bdr);

// Mark both as essential
Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 0;
ess_bdr[1-1] = 1;
ess_bdr[3-1] = 1;
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
```

### Vector Problems (H(curl) Space - ND Elements)

For electromagnetic problems, essential BCs constrain **tangential components**.

#### PEC (Perfect Electric Conductor): $n \times E = 0$

```cpp
ND_FECollection fec(order, dim);
FiniteElementSpace fes(&mesh, &fec);

Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 0;
ess_bdr[pec_attr - 1] = 1;  // PEC walls

Array<int> ess_tdof_list;
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

GridFunction E(&fes);
E = 0.0;  // Zero tangential field

a.FormLinearSystem(ess_tdof_list, E, b, A, X, B);
```

#### Port Excitation: $n \times E = f(x)$ (Modal Field)

```cpp
// Define modal field (e.g., TE10 in rectangular waveguide)
void te10_field(const Vector &x, Vector &E) {
    double a = 0.02286;  // waveguide width
    E.SetSize(3);
    E(0) = 0.0;
    E(1) = sin(M_PI * x(0) / a);  // Ey component
    E(2) = 0.0;
}

VectorFunctionCoefficient modal_field(3, te10_field);

Array<int> port_bdr(mesh.bdr_attributes.Max());
port_bdr = 0;
port_bdr[port_attr - 1] = 1;

GridFunction E(&fes);
E = 0.0;
E.ProjectBdrCoefficientTangent(modal_field, port_bdr);

Array<int> ess_bdr = port_bdr;  // Mark port as essential
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
```

### Vector Problems (H(div) Space - RT Elements)

For fluid/porous media problems, essential BCs constrain **normal components**.

#### No Penetration: $u \cdot n = 0$

```cpp
RT_FECollection fec(order, dim);
FiniteElementSpace fes(&mesh, &fec);

Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 0;
ess_bdr[wall_attr - 1] = 1;

Array<int> ess_tdof_list;
fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

GridFunction u(&fes);
u = 0.0;  // Zero normal component
```

#### Prescribed Normal Flux: $u \cdot n = g(x)$

```cpp
// Define normal flux function
void flux_function(const Vector &x, Vector &u) {
    // Normal component will be extracted by ProjectBdrCoefficientNormal
    u.SetSize(3);
    u(0) = 0.0;
    u(1) = 0.0;
    u(2) = 1.0;  // Flow in z-direction
}

VectorFunctionCoefficient flux_vec(3, flux_function);

Array<int> inlet_bdr(mesh.bdr_attributes.Max());
inlet_bdr = 0;
inlet_bdr[inlet_attr - 1] = 1;

GridFunction u(&fes);
u = 0.0;
u.ProjectBdrCoefficientNormal(flux_vec, inlet_bdr);
```

## Natural (Neumann) Boundary Conditions

Natural BCs are enforced weakly via boundary integrals added to the linear form (RHS).

### Scalar Problems: $\nabla u \cdot n = q$ on $\Gamma_N$

Weak form contribution:
$$\int_{\Gamma_N} q \, v \, ds$$

```cpp
// Heat flux: -k ∇u·n = q  =>  ∇u·n = -q/k
double flux_value = 5000.0;  // W/m^2 (heat entering domain)
ConstantCoefficient q(flux_value);

Array<int> neumann_bdr(mesh.bdr_attributes.Max());
neumann_bdr = 0;
neumann_bdr[flux_attr - 1] = 1;

LinearForm b(&fes);
b.AddBoundaryIntegrator(new BoundaryLFIntegrator(q), neumann_bdr);
b.Assemble();
```

**Sign convention:**

- Positive $q$ typically means flux **into** the domain
- Check weak form sign conventions in your specific problem

### Vector Problems: Traction on $\Gamma_N$

For elasticity: prescribed traction $t = \sigma \cdot n$

```cpp
void traction_vector(const Vector &x, Vector &t) {
    t.SetSize(3);
    t(0) = 0.0;
    t(1) = -1000.0;  // Downward force [N/m^2]
    t(2) = 0.0;
}

VectorFunctionCoefficient traction(3, traction_vector);

Array<int> traction_bdr(mesh.bdr_attributes.Max());
traction_bdr = 0;
traction_bdr[surface_attr - 1] = 1;

LinearForm b(&fes);
b.AddBoundaryIntegrator(new VectorBoundaryLFIntegrator(traction), traction_bdr);
b.Assemble();
```

## Robin (Mixed) Boundary Conditions

Robin BCs couple solution value and flux: $\alpha u + \beta \nabla u \cdot n = g$ on $\Gamma_R$

Common form: $-k \frac{\partial u}{\partial n} = h(u - u_{\infty})$ (convective heat transfer)

### Scalar Problems: Convective BC

Rearranged: $-k \frac{\partial u}{\partial n} = h \, u - h \, u_{\infty}$

Weak form:

$$\int_{\Gamma_R} h \, u \, v \, ds = \int_{\Gamma_R} h \, u_{\infty} \, v \, ds$$

- LHS term (with $u$) goes into bilinear form
- RHS term (constant) goes into linear form

```cpp
// Parameters
double h = 100.0;      // Heat transfer coefficient [W/m^2/K]
double T_inf = 300.0;  // Ambient temperature [K]

ConstantCoefficient h_coef(h);
ConstantCoefficient rhs_coef(h * T_inf);

Array<int> robin_bdr(mesh.bdr_attributes.Max());
robin_bdr = 0;
robin_bdr[convective_attr - 1] = 1;

// Bilinear form: add h*u*v on boundary
BilinearForm a(&fes);
a.AddDomainIntegrator(new DiffusionIntegrator(k));  // Volume term
a.AddBoundaryIntegrator(new MassIntegrator(h_coef), robin_bdr);  // Boundary term
a.Assemble();

// Linear form: add h*T_inf*v on boundary
LinearForm b(&fes);
b.AddBoundaryIntegrator(new BoundaryLFIntegrator(rhs_coef), robin_bdr);
b.Assemble();
```

### Impedance BC (Electromagnetics)

For time-harmonic Maxwell: $n \times (\nabla \times E) = -jk Z \, (n \times E)$

Where $Z$ is surface impedance. This is a Robin-type BC for vector fields.

```cpp
// Surface impedance BC (requires complex arithmetic or approximation)
double Z = 377.0;  // Free-space impedance [Ohms]
double k = 2.0 * M_PI * freq / c0;

// For real-valued approximation:
ConstantCoefficient impedance_coef(k * Z);
Array<int> impedance_bdr(mesh.bdr_attributes.Max());
impedance_bdr = 0;
impedance_bdr[absorbing_attr - 1] = 1;

BilinearForm a(&fes);
// ... curl-curl and mass terms ...
a.AddBoundaryIntegrator(new VectorFEMassIntegrator(impedance_coef), impedance_bdr);
a.Assemble();
```

## Periodic Boundary Conditions

Periodic BCs: $u(x + L) = u(x)$ where $L$ is the period.

### Method 1: Periodic Mesh

Use meshes with built-in periodicity (e.g., `periodic-square.mesh`).

```cpp
Mesh mesh("periodic-square.mesh");
// Periodic BCs are automatically handled in the finite element space
H1_FECollection fec(order, dim);
FiniteElementSpace fes(&mesh, &fec);
// No explicit BC marking needed
```

### Method 2: Constraint Enforcement (Manual)

For non-periodic meshes, manually enforce periodicity by identifying and constraining DOF pairs.

```cpp
// Find DOFs on opposite boundaries (attr 1 and attr 3, for example)
// This requires careful geometric matching and is problem-specific
// ... (Advanced topic, see MFEM examples with periodic constraints)
```

## Symmetry and Anti-Symmetry

### Symmetry Plane

**For scalars:** Neumann BC $\nabla u \cdot n = 0$ (natural, no flux)

```cpp
// No explicit enforcement needed; homogeneous Neumann is the default
// Just don't mark the symmetry boundary as essential or natural
```

**For vectors:** Depends on field type

- $n \times E = 0$ for electric field (use essential BC with zero tangential)
- $E \cdot n = 0$ for other fields (use essential BC with zero normal)

### Anti-Symmetry Plane

**For scalars:** Dirichlet $u = 0$ (essential)

```cpp
Array<int> ess_bdr(mesh.bdr_attributes.Max());
ess_bdr = 0;
ess_bdr[antisymmetry_attr - 1] = 1;

GridFunction u(&fes);
u = 0.0;  // Zero on anti-symmetry plane
```

## Time-Dependent Boundary Conditions

For transient problems, BCs may vary with time.

### Time-Dependent Dirichlet

```cpp
class TimeDependentDirichlet : public Coefficient {
    double t_current;
public:
    void SetTime(double t) { t_current = t; }
    
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
        Vector x;
        T.Transform(ip, x);
        // Example: sinusoidal excitation
        return sin(2.0 * M_PI * freq * t_current);
    }
};

TimeDependentDirichlet time_bc;

// In time loop:
for (int ti = 0; ti < num_steps; ti++) {
    time_bc.SetTime(t);
    u.ProjectBdrCoefficient(time_bc, ess_bdr);
    // ... solve time step ...
    t += dt;
}
```

### Time-Dependent Neumann

```cpp
class TimeDependentFlux : public Coefficient {
    double t_current;
public:
    void SetTime(double t) { t_current = t; }
    
    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip) {
        return flux_amplitude * exp(-t_current / tau);  // Decaying flux
    }
};

TimeDependentFlux flux_bc;

// In time loop:
for (int ti = 0; ti < num_steps; ti++) {
    flux_bc.SetTime(t);
    b.Update();  // Clear previous
    b.AddBoundaryIntegrator(new BoundaryLFIntegrator(flux_bc), neumann_bdr);
    b.Assemble();
    // ... solve time step ...
    t += dt;
}
```

## Advanced: Weak Enforcement of Dirichlet BCs (Nitsche's Method)

Instead of strong enforcement, add penalty term to bilinear form:

$$a(u,v) + \int_{\Gamma_D} \left( -\nabla u \cdot n \, v - u \, \nabla v \cdot n + \frac{\gamma}{h} u \, v \right) ds$$

Useful for:

- Cut-cell / immersed boundary methods
- Unfitted FEM

```cpp
// Nitsche penalty parameter
double gamma = 10.0 * order * order;  // Depends on element order
double h_mesh = mesh.GetElementSize(0);  // Representative element size

// Penalty term: (gamma/h) * u * v
ConstantCoefficient penalty_coef(gamma / h_mesh);
a.AddBoundaryIntegrator(new MassIntegrator(penalty_coef), dirichlet_bdr);

// Additional consistency terms (see Nitsche's method literature)
// ... (Advanced topic)
```

## Common Pitfalls and Debugging

### 1. Attribute Indexing

**Pitfall:** Confusing 1-indexed mesh attributes with 0-indexed arrays.

```cpp
// WRONG: using attribute value directly as index
bdr_marker[attr] = 1;  // Will access out of bounds if attr > max

// CORRECT: subtract 1
bdr_marker[attr - 1] = 1;
```

### 2. Missing FormLinearSystem Call

**Pitfall:** Forgetting to use `FormLinearSystem` which enforces essential BCs.

```cpp
// WRONG:
SparseMatrix &A = a.SpMat();
solver.SetOperator(A);
solver.Mult(b, x);

// CORRECT:
SparseMatrix A_constrained;
Vector B, X;
a.FormLinearSystem(ess_tdof_list, x, b, A_constrained, X, B);
solver.SetOperator(A_constrained);
solver.Mult(B, X);
a.RecoverFEMSolution(X, b, x);
```

### 3. Sign Convention Errors

Different PDE formulations use different sign conventions. Always check:

- Is $-\nabla^2 u$ or $+\nabla^2 u$ in your PDE?
- Is flux positive inward or outward?

Verify with simple test cases (e.g., 1D problem with known analytical solution).

### 4. BC Compatibility

For elliptic problems, pure Neumann BC requires compatibility:

$$\int_\Omega f \, dx + \int_{\Gamma_N} q \, ds = 0$$

Otherwise, problem is ill-posed (or has non-unique solution up to constant).

### 5. Forgetting to Update Time-Dependent BCs

In time loops, always update coefficients and re-project/re-assemble.

## Summary Checklist

- [ ] Identify boundary regions and check mesh attributes
- [ ] Classify each boundary as essential, natural, or Robin
- [ ] For essential: mark attributes and project boundary values
- [ ] For natural: add boundary integrators to linear form
- [ ] For Robin: add boundary integrators to both bilinear and linear forms
- [ ] Use `FormLinearSystem` to enforce essential BCs in matrix system
- [ ] Verify sign conventions match your PDE formulation
- [ ] For time-dependent: update coefficients in time loop
- [ ] Test with simple geometries and analytical solutions

## References

- MFEM examples: All `ex*.cpp` files demonstrate various BCs
- Evans, L.C., "Partial Differential Equations" (theory)
- MFEM documentation: Boundary integrators and coefficient classes
- Nitsche's method: Freund & Stenberg (1995)
