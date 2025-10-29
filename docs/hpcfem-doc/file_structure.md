# hpcfem source file structure

This document describes the refactored `src/hpcfem` layout used by the project.

Current layout (important files):

```bash
src/
└── hpcfem/
    ├── core/
    │   ├── fem_problem.hpp
    │   ├── fem_problem.cpp
    │   ├── physics_interface.hpp
    │   └── solver_interface.hpp
    ├── physics/
    │   ├── physics_electrostatics.hpp
    │   ├── physics_electrostatics.cpp
    │   ├── physics_thermal.hpp
    │   ├── physics_thermal.cpp
    │   ├── physics_joule_heating.hpp
    │   └── physics_joule_heating.cpp
    └── solvers/
        ├── solver_hypre_amg.hpp
        ├── solver_hypre_amg.cpp
```

Notes and rationale

- Files are organized by responsibility: `core` for central interfaces and orchestration, `physics` for physics modules, and `solvers` for solver implementations.
- Per project guidelines, all includes should use paths relative to `src/` (e.g. `#include "hpcfem/core/fem_problem.hpp"`).
- The refactor preserves external behavior; tests and examples were updated to use new header paths.

If you need the original file list or further re-organization (e.g. split `integrators/` or add `utils/`), we can plan a follow-up PR.
