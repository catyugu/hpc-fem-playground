# MFEM Plugin System Guide

## Overview

The MFEM Plugin System provides a modular, extensible framework for optimizing and customizing finite element method workflows. It allows you to easily swap assembly strategies and linear solvers to achieve optimal performance for your specific problem characteristics.

## Key Features

- **Drop-in replacement**: Plugins work with existing MFEM code
- **Performance instrumentation**: Automatic timing and metrics collection
- **Type-safe interfaces**: Compile-time checking of plugin compatibility
- **Header-only library**: No additional linking required
- **Extensible design**: Easy to add new strategies

## Architecture

### Plugin Base Classes

All plugins inherit from `PluginBase` which provides:

- Performance metrics tracking
- Timing utilities
- Plugin identification

```cpp
struct PerformanceMetrics {
    double setup_time_ms;      // Setup/initialization time
    double solve_time_ms;      // Solve/computation time
    double total_time_ms;      // Total time
    int iterations;            // Iterations (for solvers)
    double residual;           // Final residual
    size_t memory_bytes;       // Memory usage estimate
};
```

### Plugin Types

1. **Assembly Plugins** (`assembly_plugin.hpp`)
   - Control how bilinear forms are assembled
   - Trade-offs between memory and computation

2. **Solver Plugins** (`solver_plugin.hpp`)
   - Linear solver strategies
   - Preconditioner options

## Assembly Strategies

### Full Assembly (`AssemblyStrategy::FULL`)

Standard MFEM assembly creating a full sparse matrix.

**When to use:**

- Small to medium problems
- Multiple matrix-vector products needed
- Direct solvers
- Memory is not a constraint

**Example:**

```cpp
auto assembly = create_assembly_plugin(AssemblyStrategy::FULL);
BilinearForm a(&fes);
a.AddDomainIntegrator(new DiffusionIntegrator(coeff));
assembly->assemble(a);
```

**Performance characteristics:**

- Setup time: Moderate (one-time assembly)
- Apply time: Fast (direct SpMV)
- Memory: High (full matrix storage)

### Partial Assembly (`AssemblyStrategy::PARTIAL`)

Assembles at element level, recomputes contributions during matvec.

**When to use:**

- Large problems with memory constraints
- High-order finite elements
- Few matrix-vector products

**Example:**

```cpp
auto assembly = create_assembly_plugin(AssemblyStrategy::PARTIAL);
BilinearForm a(&fes);
a.AddDomainIntegrator(new DiffusionIntegrator(coeff));
assembly->assemble(a);
```

**Performance characteristics:**

- Setup time: Low (element matrix storage)
- Apply time: Moderate (recompute + accumulate)
- Memory: Medium (element matrices only)

### Matrix-Free (`AssemblyStrategy::MATRIX_FREE`)

No matrix storage; recomputes everything during matvec.

**When to use:**

- Very large problems
- Extreme memory constraints
- Time-dependent problems (coefficients change)

**Example:**

```cpp
auto assembly = create_assembly_plugin(AssemblyStrategy::MATRIX_FREE);
BilinearForm a(&fes);
a.AddDomainIntegrator(new DiffusionIntegrator(coeff));
assembly->assemble(a);
```

**Performance characteristics:**

- Setup time: Minimal
- Apply time: Slow (full recomputation)
- Memory: Minimal (no matrix)

### Element-by-Element (`AssemblyStrategy::ELEMENT_BY_ELEMENT`)

Custom storage of element matrices with optimized access patterns.

**When to use:**

- Custom optimization strategies
- Research on assembly algorithms
- GPU acceleration (future)

**Example:**

```cpp
auto assembly = create_assembly_plugin(AssemblyStrategy::ELEMENT_BY_ELEMENT);
BilinearForm a(&fes);
a.AddDomainIntegrator(new DiffusionIntegrator(coeff));
assembly->assemble(a);
```

## Solver Strategies

### Conjugate Gradient (`SolverType::CG`)

Krylov method for symmetric positive definite systems.

**When to use:**

- Poisson/diffusion problems
- Elasticity problems
- Most H1 problems

**Preconditioners:**

- `NONE`: Unpreconditioned CG (slow convergence)
- `JACOBI`: Diagonal scaling (cheap, modest improvement)
- `GAUSS_SEIDEL`: Forward/backward sweep (better than Jacobi)
- `AMG`: Algebraic multigrid (best for large problems)

**Example:**

```cpp
auto solver = create_solver_plugin(SolverType::CG, PreconditionerType::JACOBI);
solver->set_tolerance(1e-9, 1e-12);
solver->set_max_iterations(1000);
solver->setup(A);
solver->solve(B, X);
```

### GMRES (`SolverType::GMRES`)

Generalized minimal residual for non-symmetric systems.

**When to use:**

- Convection-diffusion problems
- Non-symmetric formulations
- Saddle-point problems

**Options:**

```cpp
auto solver = create_solver_plugin(SolverType::GMRES, PreconditionerType::JACOBI);
// GMRES specific: Krylov subspace dimension (default: 50)
solver->set_max_iterations(1000);
```

### AMG-CG (`SolverType::AMG_CG`)

CG with algebraic multigrid preconditioner.

**When to use:**

- Large 3D problems
- Complex geometries
- High mesh refinement

**Example:**

```cpp
auto solver = create_solver_plugin(SolverType::AMG_CG, PreconditionerType::AMG);
solver->set_tolerance(1e-9, 1e-12);
solver->setup(A);
solver->solve(B, X);
```

### Direct Solver (`SolverType::DIRECT`)

UMFPACK direct solver (requires SuiteSparse).

**When to use:**

- Small problems (< 100k DOFs)
- Multiple RHS vectors
- Need exact solution

**Example:**

```cpp
#ifdef MFEM_USE_SUITESPARSE
auto solver = create_solver_plugin(SolverType::DIRECT, PreconditionerType::NONE);
solver->setup(A);
solver->solve(B, X);
#endif
```

## Complete Workflow Example

### Basic Usage

```cpp
#include "mfem.hpp"
#include "plugins/mfem_plugins.hpp"

using namespace mfem;
using namespace mfem_plugins;

int main() {
    // Setup problem (mesh, FE space, etc.)
    Mesh mesh("mesh.msh");
    H1_FECollection fec(2, mesh.Dimension());
    FiniteElementSpace fes(&mesh, &fec);
    
    // Choose assembly strategy
    auto assembly = create_assembly_plugin(AssemblyStrategy::FULL);
    
    // Assemble bilinear form
    BilinearForm a(&fes);
    ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    assembly->assemble(a);
    
    // Setup RHS and BCs
    LinearForm b(&fes);
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();
    
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1;
    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    
    GridFunction x(&fes);
    x = 0.0;
    
    // Form system
    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
    
    // Choose solver strategy
    auto solver = create_solver_plugin(SolverType::CG, PreconditionerType::JACOBI);
    solver->set_tolerance(1e-9, 1e-12);
    solver->set_max_iterations(1000);
    
    // Solve
    solver->setup(A);
    solver->solve(B, X);
    
    // Recover solution
    a.RecoverFEMSolution(X, b, x);
    
    // Get metrics
    std::cout << "Assembly time: " << assembly->metrics().setup_time_ms << " ms\n";
    std::cout << "Solver time: " << solver->metrics().solve_time_ms << " ms\n";
    std::cout << "Iterations: " << solver->metrics().iterations << "\n";
    
    return 0;
}
```

### Runtime Configuration

```cpp
// Parse command-line options
AssemblyStrategy asm_strategy = AssemblyStrategy::FULL;
SolverType solver_type = SolverType::CG;
PreconditionerType precond = PreconditionerType::JACOBI;

for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "--assembly") == 0) {
        std::string type = argv[++i];
        if (type == "full") asm_strategy = AssemblyStrategy::FULL;
        else if (type == "partial") asm_strategy = AssemblyStrategy::PARTIAL;
        else if (type == "matfree") asm_strategy = AssemblyStrategy::MATRIX_FREE;
    }
    // ... similar for solver and preconditioner
}

// Use selected strategies
auto assembly = create_assembly_plugin(asm_strategy);
auto solver = create_solver_plugin(solver_type, precond);
```

## Performance Tuning Guidelines

### Problem Size-Based Recommendations

**Small problems (< 10k DOFs):**

```cpp
auto assembly = create_assembly_plugin(AssemblyStrategy::FULL);
auto solver = create_solver_plugin(SolverType::DIRECT, PreconditionerType::NONE);
```

**Medium problems (10k - 100k DOFs):**

```cpp
auto assembly = create_assembly_plugin(AssemblyStrategy::FULL);
auto solver = create_solver_plugin(SolverType::CG, PreconditionerType::GAUSS_SEIDEL);
```

**Large problems (100k - 1M DOFs):**

```cpp
auto assembly = create_assembly_plugin(AssemblyStrategy::PARTIAL);
auto solver = create_solver_plugin(SolverType::AMG_CG, PreconditionerType::AMG);
```

**Very large problems (> 1M DOFs):**

```cpp
auto assembly = create_assembly_plugin(AssemblyStrategy::MATRIX_FREE);
auto solver = create_solver_plugin(SolverType::AMG_CG, PreconditionerType::AMG);
```

### Geometry-Based Recommendations

**Simple geometries (cubes, cylinders):**

- Full assembly + CG with Jacobi
- Direct solver for small problems

**Complex geometries:**

- Full assembly for accuracy
- AMG preconditioner for robustness

**Extreme aspect ratios:**

- AMG preconditioner essential
- May need custom coarsening strategies

### Time-Dependent Problems

For problems with changing coefficients:

```cpp
// Matrix-free is efficient if coefficients change
auto assembly = create_assembly_plugin(AssemblyStrategy::MATRIX_FREE);

// Time loop
for (double t = 0; t < t_final; t += dt) {
    // Update coefficients
    coeff.SetTime(t);
    
    // Reassemble (fast with matrix-free)
    assembly->assemble(a);
    
    // Solve
    solver->solve(B, X);
}
```

## Benchmarking

Use the provided benchmarks to find optimal settings:

### Assembly Benchmark

```bash
./bench_assembly testdata/testmesh_cube.mesh 2 2
```

### Solver Benchmark

```bash
./bench_solver testdata/testmesh_cube.mesh 2 2
```

### Combined Workflow

```bash
./bench_combined testdata/testmesh_cube.mesh 2 2
```

## Extending the Plugin System

### Adding Custom Assembly Strategy

Create new enum value:

```cpp
enum class AssemblyStrategy {
    FULL,
    PARTIAL,
    MATRIX_FREE,
    ELEMENT_BY_ELEMENT,
    MY_CUSTOM_STRATEGY  // Add here
};
```

Implement plugin class:

```cpp
class MyCustomAssemblyPlugin : public AssemblyPlugin {
public:
    void assemble(BilinearForm& a) override {
        // Custom assembly logic
    }
    
    void mult(const Vector& x, Vector& y) const override {
        // Custom matvec
    }
    
    size_t memory_usage() const override {
        // Estimate memory
        return 0;
    }
};
```

Update factory function:

```cpp
std::unique_ptr<AssemblyPlugin> create_assembly_plugin(AssemblyStrategy strategy) {
    switch (strategy) {
        // ... existing cases
        case AssemblyStrategy::MY_CUSTOM_STRATEGY:
            return std::make_unique<MyCustomAssemblyPlugin>();
    }
}
```

### Adding Custom Solver

Follow similar pattern in `solver_plugin.hpp`.

## Best Practices

1. **Profile before optimizing**: Use benchmarks to identify bottlenecks
2. **Start simple**: Begin with full assembly + preconditioned CG
3. **Measure everything**: Use plugin metrics for quantitative decisions
4. **Test convergence**: Verify solution accuracy with different strategies
5. **Document choices**: Comment why specific plugins were selected

## Common Issues

### Solver doesn't converge

- Try stronger preconditioner (GAUSS_SEIDEL or AMG)
- Increase max iterations
- Check problem conditioning
- Verify boundary conditions

### Out of memory

- Switch to partial or matrix-free assembly
- Reduce mesh refinement
- Use iterative solver instead of direct

### Slow performance

- Use AMG for large problems
- Consider matrix-free for many matvecs
- Check if problem is well-conditioned
- Profile to find bottleneck (assembly vs solve)

## References

- MFEM documentation: <https://mfem.org/>
- See `benchmark/` directory for performance comparisons
- See `example/ex5_plugin_demo.cpp` for complete working example
- See individual plugin headers for API details

## Version

Plugin system version: 1.0.0
Compatible with: MFEM 4.x
