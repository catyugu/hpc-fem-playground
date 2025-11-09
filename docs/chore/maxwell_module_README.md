# Maxwell Time-Domain Physics Module

## Overview

This module implements a time-domain solver for the coupled first-order Maxwell equations:

```
ε ∂E/∂t = ∇×(μ⁻¹B) - σE - J
∂B/∂t = -∇×E
```

## Features

- **Mixed Formulation**: H(curl) discretization for E-field, H(div) for B-field
- **Energy Conservation**: Symplectic time integration preserves energy in lossless case
- **Material Support**: Dielectric, magnetic, and conductive materials
- **Source Terms**: Time-dependent volumetric current densities J(x,t)
- **Boundary Conditions**:
  - Natural (zero tangential current)
  - Dirichlet (prescribed ∂E/∂t)
  - Sommerfeld ABC (first-order absorbing)
  - TODO: PML, Port BCs
- **MPI Parallel**: Full parallel support via MFEM/HYPRE

## Quick Start

```cpp
#include "hpcfem/physics/physics_maxwell_timedomain.hpp"

// Define material properties
double epsilon(const Vector& x) { return 8.854e-12; }  // Permittivity
double muInv(const Vector& x) { return 1.0/(4e-7*M_PI); }  // μ⁻¹
double sigma(const Vector& x) { return 0.0; }  // Lossless

// Define current source (optional)
void currentSource(const Vector& x, double t, Vector& j) {
    j.SetSize(3);
    j = 0.0;
    // Define your source here
}

// Create solver
PhysicsMaxwellTimeDomain solver(
    pmesh,           // Parallel mesh
    2,               // Polynomial order
    epsilon,         // ε(x)
    muInv,           // μ⁻¹(x)
    sigma,           // σ(x)
    currentSource,   // J(x,t)
    abcMarkers,      // ABC boundaries
    dbcMarkers,      // Dirichlet boundaries
    nullptr          // ∂E/∂t BC function
);

// Set initial conditions
VectorFunctionCoefficient e0(...);
solver.setInitialEField(e0);

// Time stepping (use MFEM ODE solvers)
RK4Solver ode_solver;
ode_solver.Init(solver);
double t = 0.0;
double dt = solver.getMaximumTimeStep() * 0.5;
for (int step = 0; step < numSteps; step++) {
    ode_solver.Step(B, t, dt);
    t += dt;
    
    // Monitor energy
    double energy = solver.getEnergy();
    cout << "t = " << t << ", E = " << energy << endl;
}
```

## Architecture

The module is split into 4 components:

1. **PhysicsMaxwellTimeDomain** - Main time-dependent operator
2. **MaxwellMaterialProperties** - Material coefficient management
3. **MaxwellSourceTerms** - Current density sources
4. **MaxwellBoundaryConditions** - BC handling

## Current Status

✅ **Complete:**
- Architecture and interfaces
- Material, source, and BC management
- Matrix assembly
- MPI parallel support
- Test framework

⏳ **In Progress:**
- Operator evaluation (`Mult()`)
- Energy computation
- Implicit solver setup
- Test implementations

🔮 **Future:**
- PML boundary conditions
- Port boundary conditions
- Higher-order time integration (orders 2-4)
- Frequency-domain post-processing

## Testing

```bash
# Build and run tests
cd cmake-build-debug
make test_physics_maxwell_timedomain
mpirun -np 4 ./tests/test_physics_maxwell_timedomain
```

Current tests (framework ready, awaiting implementation):
1. Energy conservation in lossless cavity
2. Current source radiation
3. Lossy material energy decay
4. Absorbing boundary conditions
5. MPI parallel correctness

## Examples

TODO: Create examples
- `ex7_maxwell_cavity.cpp` - Cavity resonance
- `ex8_maxwell_waveguide_pulse.cpp` - Pulse propagation  
- `ex9_maxwell_lossy.cpp` - Lossy materials

## References

- MFEM maxwell miniapp: `share/mfem/miniapps/electromagnetics/`
- Maxwell's equations: [Wikipedia](https://en.wikipedia.org/wiki/Maxwell%27s_equations)
- MFEM finite elements: [MFEM documentation](https://mfem.org)

## Contributing

Follow TDD workflow:
1. Write/update tests
2. Implement functionality
3. Verify tests pass
4. Refactor and document

Adhere to project guidelines in `.github/instructions/RULES.instructions.md`.

## License

Same as parent project (see top-level LICENSE).
