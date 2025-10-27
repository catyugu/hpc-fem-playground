# HPC FEM Playground - MFEM Learning Framework

## Overview

This repository provides a comprehensive learning and development framework for MFEM (Modular Finite Element Methods library), focusing on:

1. **Comprehensive Documentation** - Detailed guides for different problem types
2. **Working Examples** - Ready-to-run examples demonstrating various PDE problems
3. **Extensible Plugin System** - Modular components for optimization and algorithm research

## Documentation (`docs/mfem-design-concept/`)

Comprehensive startup guides covering:

### 1. Static Problems (`basic_usage_scalar_static.md`)

- Steady-state scalar diffusion/heat problems
- Pure Dirichlet and mixed boundary conditions
- Step-by-step implementation patterns

### 2. Transient Problems (`transient_time_dependent.md`)

- Time-dependent parabolic (heat/diffusion) and hyperbolic (wave) equations
- ODE time integration framework (explicit and implicit methods)
- Mass and stiffness matrix assembly
- Time-dependent boundary conditions
- Stability analysis and time step selection

### 3. Frequency-Domain Problems (`frequency_domain_problems.md`)

- Time-harmonic (frequency-domain) PDEs
- Helmholtz equation for acoustics
- Vector Helmholtz / Maxwell equations for electromagnetics
- Waveguide and cavity problems
- Port excitation and modal analysis

### 4. Eigenvalue Problems (`eigenvalue_problems.md`)

- Generalized eigenvalue problems (K u = λ M u)
- Laplacian eigenvalues (vibration modes, quantum states)
- Maxwell eigenvalues (cavity modes)
- LOBPCG eigensolver with AMG preconditioning
- Validation against analytical solutions

### 5. Boundary Conditions (`boundary_conditions_guide.md`)

- Essential (Dirichlet) BCs for scalar and vector fields
- Natural (Neumann) BCs
- Robin (mixed) BCs
- Periodic BCs
- Time-dependent BCs
- Boundary attribute marking and management

## Examples (`example/`)

### ex0_cylinder_dirichlet.cpp

- **Problem**: Steady-state heat conduction in cylinder
- **Type**: Static, pure Dirichlet BC
- **Features**: Multiple temperature boundaries

### ex1_cube_mixed_bc.cpp

- **Problem**: Steady-state heat with mixed BCs
- **Type**: Static, Dirichlet + Robin + Neumann
- **Features**: Convective heat transfer, prescribed flux

### ex2_rect_waveguide.cpp

- **Problem**: Rectangular waveguide electromagnetics
- **Type**: Frequency-domain, H(curl) space
- **Features**: TE10 modal excitation, PEC walls

### ex3_transient_heat.cpp ✓ NEW

- **Problem**: Time-dependent heat diffusion
- **Type**: Transient, parabolic PDE
- **Features**:
  - Backward Euler, Forward Euler, and RK4 time integrators
  - Gaussian initial condition
  - Dirichlet BC
  - ParaView time series output

### ex4_eigenvalue_laplacian.cpp ✓ NEW

- **Problem**: Laplacian eigenvalue problem (-∇²u = λu)
- **Type**: Eigenvalue problem
- **Features**:
  - LOBPCG solver with BoomerAMG preconditioner
  - Homogeneous Dirichlet BC
  - Validation against analytical solutions (unit cube)
  - Multiple eigenmode extraction and visualization

## Building and Running

```bash
# Configure
cmake -B cmake-build-debug -DCMAKE_BUILD_TYPE=Debug

# Build all examples
cmake --build cmake-build-debug -j4

# Run examples
./cmake-build-debug/example/ex3_transient_heat testdata/testmesh_cube.mesh results/ex3 1 0.1 0.01 2
./cmake-build-debug/example/ex4_eigenvalue_laplacian testdata/testmesh_cube.mesh results/ex4 1 5

# Visualize in ParaView
paraview results/ex3/ex3.pvd
paraview results/ex4/ex4.pvd
```

## Test Meshes (`testdata/`)

Available meshes for testing:

- `testmesh_bar.mesh` - 3D bar/beam geometry
- `testmesh_cube.mesh` - 3D cube with 6 boundary attributes
- `testmesh_cylinder.mesh` - 3D cylinder geometry

All meshes include properly marked boundary attributes for different BC types.

## Plugin System (In Development - `src/`)

Modular, extensible components for replacing MFEM pipeline stages:

### Planned Components

1. **Assembly Plugins**
   - Element-by-element assembly
   - Matrix-free operators
   - Partial assembly strategies
   - Performance benchmarking hooks

2. **Solver Plugins**
   - Custom preconditioners
   - Domain decomposition methods
   - Multigrid variants
   - Iterative solver wrappers with metrics

3. **Time Integrator Plugins**
   - Custom ODE solvers
   - Adaptive time stepping
   - Multi-rate methods

## Usage Workflow

### For Learning MFEM

1. Read the relevant guide in `docs/mfem-design-concept/`
2. Examine the corresponding example in `example/`
3. Modify and experiment with parameters
4. Visualize results in ParaView

### For Research/Development

1. Start with existing examples as templates
2. Implement custom algorithms using plugin interfaces (when available)
3. Benchmark performance against standard methods
4. Iterate and optimize

## Key MFEM Concepts Demonstrated

- **Finite Element Spaces**: H1, H(curl), H(div)
- **Bilinear Forms**: Stiffness and mass matrices
- **Linear Forms**: Load vectors and boundary integrals
- **Grid Functions**: Solution representation
- **Integrators**: Domain and boundary term assembly
- **Solvers**: CG, GMRES, direct solvers (UMFPACK, MUMPS)
- **Preconditioners**: AMG, AMS, smoothers
- **Time Integration**: Explicit and implicit ODE solvers
- **Eigensolvers**: LOBPCG for large-scale problems

## Dependencies

- MFEM (included in `share/mfem/`)
- HYPRE (for parallel solvers and AMG)
- MPI (optional, for parallel execution)
- SuiteSparse/UMFPACK (optional, for direct solvers)
- ParaView (for visualization)

## Next Steps

1. ✅ Documentation - Complete
2. ✅ Basic examples - Complete
3. ✅ Transient example - Complete
4. ✅ Eigenvalue example - Complete
5. 🚧 Plugin architecture - In Progress
6. ⏳ Assembly plugins - Planned
7. ⏳ Solver plugins - Planned
8. ⏳ Performance benchmarks - Planned

## References

- [MFEM Website](https://mfem.org)
- [MFEM Examples](https://mfem.org/examples/)
- [MFEM Documentation](https://docs.mfem.org)
- HYPRE documentation
- Finite Element Method textbooks

## Contributing

This is a learning and research framework. Feel free to:

- Add new problem types and examples
- Implement optimization plugins
- Improve documentation
- Report issues or suggest improvements

## License

Follow MFEM's licensing (BSD-3-Clause for MFEM core).
