# Library Structure and Organization

## Overview

This document explains the library architecture and how different components are linked together in the HPC FEM Playground project.

## Library Hierarchy

```bash
mfem (external library)
  ↓
mfem_plugins (header-only, src/)
  ↓
HPC_FEM_PLAYGROUND_LIB (INTERFACE library)
  ↓
examples/ and benchmark/ executables
```

## Component Details

### 1. MFEM Library

- **Location**: Built in `cmake-build-*/mfem_build/`
- **Type**: Static library
- **Purpose**: Core finite element library providing mesh, finite element spaces, operators, solvers
- **Dependencies**: HYPRE, SuiteSparse, METIS, MPI (optional)

### 2. mfem_plugins

- **Location**: `src/plugins/`
- **Type**: INTERFACE library (header-only)
- **Purpose**: Performance optimization plugins for assembly strategies and linear solvers
- **Headers**:
  - `plugin_base.hpp` - Base classes, performance metrics, timing utilities
  - `assembly_plugin.hpp` - Assembly strategy plugins (Full, Partial, Matrix-Free, Element-by-Element)
  - `solver_plugin.hpp` - Linear solver plugins (CG, GMRES, Direct) with preconditioners
  - `mfem_plugins.hpp` - Main header including all plugins
- **Dependencies**: None directly (will link through HPC_FEM_PLAYGROUND_LIB)

### 3. HPC_FEM_PLAYGROUND_LIB

- **Location**: Defined in root `CMakeLists.txt`
- **Type**: INTERFACE library
- **Purpose**: Unified library target combining MFEM and plugins
- **Links**: `mfem` + `mfem_plugins`
- **Includes**:
  - MFEM headers
  - Plugin headers (`${CMAKE_SOURCE_DIR}/src`)

### 4. Examples

- **Location**: `example/`
- **Executables**:
  - `ex0_cylinder_dirichlet` - Static Poisson problem
  - `ex1_cube_mixed_bc` - Mixed boundary conditions
  - `ex2_rect_waveguide` - Frequency-domain Maxwell
  - `ex3_transient_heat` - Time-dependent heat equation
  - `ex4_eigenvalue_laplacian` - Eigenvalue problem
  - `ex5_plugin_demo` - Plugin system demonstration
- **Links**: `HPC_FEM_PLAYGROUND_LIB` only

### 5. Benchmarks

- **Location**: `benchmark/`
- **Executables**:
  - `bench_assembly` - Compare assembly strategies
  - `bench_solver` - Compare solver strategies
  - `bench_combined` - Combined assembly + solver workflows
- **Links**: `HPC_FEM_PLAYGROUND_LIB` only

## Build Order

The CMake configuration ensures proper build order:

1. **MFEM** is built first (external library)
2. **mfem_plugins** INTERFACE library is created (just header references)
3. **HPC_FEM_PLAYGROUND_LIB** is created, linking both mfem and mfem_plugins
4. **Examples** and **benchmarks** are built, each linking to HPC_FEM_PLAYGROUND_LIB

## Why This Structure?

### Benefits

1. **Simplicity**: Examples and benchmarks only need to link to one target
2. **Consistency**: All executables get the same MFEM + plugins configuration
3. **Maintainability**: Changes to plugin or MFEM configuration happen in one place
4. **Modularity**: Plugins are separate from MFEM but integrated at the library level

### Design Decisions

- **Header-only plugins**: No separate compilation needed, easier to template
- **INTERFACE library for HPC_FEM_PLAYGROUND_LIB**: Transitive dependencies automatically propagate
- **Single link target**: Reduces CMake complexity in example/benchmark definitions

## Usage Pattern

### For New Examples

```cmake
add_executable(my_new_example my_new_example.cpp)
target_link_libraries(my_new_example PRIVATE HPC_FEM_PLAYGROUND_LIB)
target_include_directories(my_new_example PRIVATE 
    ${HPC_FEM_PLAYGROUND_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/src
)
```

### In Code

```cpp
#include "mfem.hpp"              // MFEM library
#include "plugins/mfem_plugins.hpp"  // Plugin system

using namespace mfem;
using namespace mfem_plugins;

int main() {
    // Use MFEM normally
    Mesh mesh("mesh.msh");
    H1_FECollection fec(2, mesh.Dimension());
    FiniteElementSpace fes(&mesh, &fec);
    
    // Use plugins for optimization
    auto assembly = create_assembly_plugin(AssemblyStrategy::FULL);
    auto solver = create_solver_plugin(SolverType::CG, PreconditionerType::JACOBI);
    
    // ... rest of code
}
```

## Directory Structure

```bash
hpc-fem-playground/
├── CMakeLists.txt                  # Main config, defines HPC_FEM_PLAYGROUND_LIB
├── src/
│   ├── CMakeLists.txt             # Defines mfem_plugins INTERFACE library
│   └── plugins/
│       ├── plugin_base.hpp         # Base classes
│       ├── assembly_plugin.hpp     # Assembly plugins
│       ├── solver_plugin.hpp       # Solver plugins
│       └── mfem_plugins.hpp        # Main plugin header
├── example/
│   ├── CMakeLists.txt             # Links executables to HPC_FEM_PLAYGROUND_LIB
│   ├── ex0_cylinder_dirichlet.cpp
│   ├── ex1_cube_mixed_bc.cpp
│   ├── ex2_rect_waveguide.cpp
│   ├── ex3_transient_heat.cpp
│   ├── ex4_eigenvalue_laplacian.cpp
│   └── ex5_plugin_demo.cpp
└── benchmark/
    ├── CMakeLists.txt             # Links executables to HPC_FEM_PLAYGROUND_LIB
    ├── bench_assembly.cpp
    ├── bench_solver.cpp
    └── bench_combined.cpp
```

## Troubleshooting

### Link Order Issues

If you see undefined references, check that `HPC_FEM_PLAYGROUND_LIB` is defined **after** `mfem_plugins`:

```cmake
add_subdirectory(src)  # Creates mfem_plugins first
add_library(HPC_FEM_PLAYGROUND_LIB INTERFACE)
target_link_libraries(HPC_FEM_PLAYGROUND_LIB INTERFACE mfem mfem_plugins)
```

### Include Path Issues

Make sure both MFEM and src directories are in include paths:

```cmake
target_include_directories(my_target PRIVATE 
    ${HPC_FEM_PLAYGROUND_INCLUDE_DIR}  # MFEM headers
    ${CMAKE_SOURCE_DIR}/src             # Plugin headers
)
```

### Plugin Availability

Plugins are header-only, so they're available as long as:

1. `${CMAKE_SOURCE_DIR}/src` is in include path
2. You `#include "plugins/mfem_plugins.hpp"`

## Performance Considerations

- **Header-only plugins**: No runtime overhead for library calls
- **Inline functions**: Compiler can optimize across boundaries
- **Template-based**: Specialization happens at compile time
- **Zero abstraction cost**: Plugin interface compiles to direct function calls

## Extension Points

To add new functionality:

1. **New assembly strategy**: Add to `assembly_plugin.hpp`, update factory function
2. **New solver**: Add to `solver_plugin.hpp`, update factory function
3. **New plugin category**: Create new header in `src/plugins/`, include in `mfem_plugins.hpp`
4. **New example**: Add .cpp to `example/`, CMake will auto-detect and build

## Version Information

- Plugin system version: 1.0.0
- MFEM version: 4.8.0
- CMake minimum: 3.16
- C++ standard: C++11 (compatible with C++14/17/20)
