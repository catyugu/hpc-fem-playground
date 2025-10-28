---
mode: agent
---

# HPC-FEM Development Guidelines and Roadmap

## **Phase 0: Project Setup and "Hello, FEM!"**

**Objective:** Establish the TDD environment, build system (CMake), and CI configuration. Ensure basic MFEM/HYPRE integration is successful and all guidelines are configured.

### **Step 0.1: Project Scaffolding and CI**

* **Planning:**
    * Create the repository structure.
    * Initialize `git` and create `dev` and `master` branches. All work will be on `dev` (Guideline #9).
    * Create a root `CMakeLists.txt` to find MFEM, HYPRE, and MPI.
    * Set up a basic CI (e.g., GitHub Actions) to run tests on every commit to `dev`.
* **Test (`tests/test_project_structure.cpp`):**
    * `test_directories_exist`: Write a test (can be a shell script run by `ctest`) that verifies the existence of `src/`, `include/hpcfem/`, `tests/`, `examples/`, `benchmarks/`, and `docs/hpcfem-doc/`.
    * `test_cmake_finds_mfem`: Write a minimal `tests/CMakeLists.txt` that configures a test. The test `test_mfem_linkage` will compile and run a C++ file that `#include "mfem.hpp"` and returns `EXIT_SUCCESS`.
* **Implementation:**
    * Create the directory structure.
    * Create `docs/hpcfem-doc/naming_registry.md`. Add `hpcfem` as the root namespace.
    * Write the root and `tests/` `CMakeLists.txt`.
* **Test:** Run `ctest`. All tests must pass.
* **Refactor/Doc:**
    * Add a `README.md` with build instructions.
    * Add `.gitignore`.
    * Commit to `dev`: `feat: Initial project setup and TDD environment`.

### **Step 0.2: Baseline MFEM Example (Replication)**

* **Planning:**
    * Replicate a standard MFEM example (e.g., Example 1) *natively* within our test framework to create a "known good" baseline. This validates the parallel build.
* **Test (`tests/test_mfem_baseline.cpp`):**
    * `test_mfem_ex1_parallel_solves_correctly`: Adapt MFEM Example 1 (serial or parallel) into a test. The test will run the solver on a simple mesh and assert that the final computed error norm is below a `constexpr` tolerance (e.g., `EXPECTED_ERROR_NORM < 1e-8`). This test *must* be runnable with MPI (`mpirun -np 4 ...`).
* **Implementation:**
    * Copy and adapt the chosen MFEM example into the test harness.
    * Ensure all constants are `constexpr` (Guideline #1, #15).
    * Add all new file names and placeholder classes/enums to `docs/hpcfem-doc/naming_registry.md` (Guideline #5).
* **Test:** Run `ctest` (both serially and with `mpirun`).
* **Refactor/Doc:**
    * Add Doxygen comments to the test file explaining its purpose (Guideline #6).
    * Commit to `dev`: `test: Replicate MFEM Ex1 as baseline test`.

---

## **Phase 1: Baseline & Replication (Joule Heating Problem)**

**Objective:** Establish a correct, validated multiphysics starting point for comparison, as per document Phase 1. We will model Joule Heating (Coupled Electrostatics and Thermal).

### **Step 1.1: Standalone Electrostatics Solver**

* **Planning:**
    * Create a standalone "physics" solver for the electrostatic problem ($\nabla \cdot (\sigma \nabla u) = f$).
    * We will use the Method of Manufactured Solutions (MMS) to verify correctness.
* **Test (`tests/test_physics_electrostatics.cpp`):**
    * `test_electrostatics_manufactured_solution`:
        1.  Define a 2D analytical solution $u(x,y) = \sin(\pi x) \cos(\pi y)$.
        2.  Compute the corresponding forcing term $f = -\nabla^2 u$.
        3.  Implement a test function that builds the FEM system ($A u = b$) for this problem using MFEM's `BilinearForm` and `LinearForm`.
        4.  Solve the system using MFEM's built-in `HypreBoomerAMG` preconditioner.
        5.  Compute the $L_2$ error norm $\|u_{fem} - u_{exact}\|$.
        6.  Assert that the error is below a `constexpr` tolerance `MMS_TOLERANCE`.
* **Implementation:**
    * Create `examples/joule_heating/solve_electrostatics.cpp`. This file will contain the logic being tested.
    * The test will call functions from this example.
    * Ensure no magic numbers (Guideline #15).
* **Test:** Run `ctest`.
* **Refactor/Doc:**
    * Clean up the example file. Add Doxygen comments.
    * Commit to `dev`: `feat: Add standalone electrostatics solver with MMS test`.

### **Step 1.2: Standalone Thermal Solver**

* **Planning:**
    * Implement a standalone solver for the heat equation ($\nabla \cdot (k \nabla T) = Q$) using MMS.
* **Test (`tests/test_physics_thermal.cpp`):**
    * `test_thermal_manufactured_solution`:
        1.  Identical to Step 1.1, but for the heat equation.
        2.  Define $T(x,y) = \cosh(x) \sinh(y)$.
        3.  Compute $Q = -\nabla^2 T$.
        4.  Build and solve the FEM system.
        5.  Assert $L_2$ error norm $\|T_{fem} - T_{exact}\| < \text{MMS_TOLERANCE}$.
* **Implementation:**
    * Create `examples/joule_heating/solve_thermal.cpp`.
* **Test:** Run `ctest`.
* **Refactor/Doc:**
    * Doxygen comments for the new example.
    * Commit to `dev`: `feat: Add standalone thermal solver with MMS test`.

### **Step 1.3: Coupled Joule Heating (One-Way)**

* **Planning:**
    * Implement a one-way coupled simulation.
    * Solve electrostatics for $u$.
    * Compute the heat source $Q = \sigma |\nabla u|^2$.
    * Solve the thermal problem for $T$ using this $Q$.
* **Test (`tests/test_physics_coupling.cpp`):**
    * `test_joule_heating_one_way_coupling_sanity_check`:
        1.  Create a simple 2D rectangular mesh.
        2.  Apply Dirichlet BCs (e.g., 1V and 0V) to the electrostatic problem.
        3.  Run the coupled simulation.
        4.  Assert that the computed `max_temperature` is greater than the ambient boundary condition temperature (a basic sanity check).
        5.  Assert that the `max_temperature` is below a reasonable upper bound (e.g., `MAX_EXPECTED_TEMP`).
* **Implementation:**
    * Create `examples/joule_heating/solve_coupled.cpp`. This script will perform the two-step solve.
* **Test:** Run `ctest`.
* **Refactor/Doc:**
    * This script is the "baseline". Add extensive comments.
    * Merge `dev` to `master`: `milestone: Phase 1 complete. Baseline Joule heating solver replicated.`.

---

## **Phase 2: Encapsulation & Abstraction**

**Objective:** Refactor the Phase 1 scripts into a modular, testable, and reusable library, as per document Phase 2 and Guideline #12 (Testability).

### **Step 2.1: Abstract Solver Interface**

* **Planning:**
    * Define a pure abstract base class for all linear solvers.
* **Test (`tests/test_solver_interface.cpp`):**
    * `test_solver_interface_compiles`: Write a test that includes the new header and defines a mock solver class inheriting from it. This is a compile-only test.
* **Implementation:**
    * Create `include/hpcfem/solver_interface.hpp`.
    * Define `namespace hpcfem { ... }` (Guideline #2).
    * Define `class SolverInterface` (Guideline #16).
    * Define `virtual void solve(const mfem::HypreParMatrix& A, const mfem::Vector& b, mfem::Vector& x) = 0;`.
    * Register `SolverInterface` in `docs/hpcfem-doc/naming_registry.md`.
* **Test:** Run `ctest`.
* **Refactor/Doc:**
    * Add Doxygen comments for the class and its methods (Guideline #6).
    * Commit to `dev`: `refactor: Introduce abstract SolverInterface`.

### **Step 2.2: Concrete HYPRE AMG Solver**

* **Planning:**
    * Wrap the `HypreBoomerAMG` logic from Phase 1 into a concrete class.
* **Test (`tests/test_solver_hypre_amg.cpp`):**
    * `test_hypre_amg_solver_solves_poisson`:
        1.  Reuse the Poisson problem setup from `test_electrostatics_manufactured_solution`.
        2.  Instantiate `hpcfem::HypreAmgSolver` (as a `SolverInterface*`).
        3.  Call `solver->solve(A, b, x)`.
        4.  Assert the $L_2$ error norm is below `MMS_TOLERANCE`.
* **Implementation:**
    * Create `include/hpcfem/solver_hypre_amg.hpp` and `src/solver_hypre_amg.cpp`.
    * `class HypreAmgSolver : public SolverInterface` (Guideline #4, no nesting).
    * Implement the `solve` method by encapsulating the `HypreBoomerAMG` setup.
    * Register new class and files in `docs/hpcfem-doc/`.
* **Test:** Run `ctest`.
* **Refactor/Doc:**
    * Ensure file length is < 500 lines (Guideline #14).
    * Commit to `dev`: `feat: Implement HypreAmgSolver wrapper`.

### **Step 2.3: Abstract Physics and Problem Classes**

* **Planning:**
    * Abstract the problem setup.
    * `Physics`: Responsible for assembling matrices/vectors.
    * `Problem`: The top-level class that *owns* a `Mesh`, `Physics`, and `Solver`.
* **Test (`tests/test_problem_abstraction.cpp`):**
    * `test_poisson_problem_solves`:
        1.  Create `hpcfem::ElectrostaticsPhysics` class.
        2.  Create `hpcfem::FemProblem` class.
        3.  The test will instantiate them:
            * `auto mesh = ...;`
            * `auto physics = hpcfem::ElectrostaticsPhysics(mesh, ...);`
            * `auto solver = hpcfem::HypreAmgSolver();`
            * `auto problem = hpcfem::FemProblem(mesh, physics, solver);`
        4.  `problem.assemble();`
        5.  `problem.solve();`
        6.  Assert the $L_2$ error norm (retrieved from `problem.getSolution()`) is below `MMS_TOLERANCE`.
* **Implementation:**
    * `include/hpcfem/physics_interface.hpp` (with `virtual void assemble(mfem::HypreParMatrix& A, mfem::Vector& b) = 0;`).
    * `include/hpcfem/physics_electrostatics.hpp`, `src/physics_electrostatics.cpp`.
    * `include/hpcfem/problem.hpp`, `src/problem.cpp`. `FemProblem` class holds `SolverInterface*`, `PhysicsInterface*`, etc.
    * Register all new classes in `docs/hpcfem-doc/`.
* **Test:** Run `ctest`.
* **Refactor:**
    * Rewrite the `examples/joule_heating/` to use this new abstract framework. The `solve_coupled.cpp` should now be much cleaner, coordinating `FemProblem` instances.
    * Ensure all examples compile and run (Guideline #7).
    * Commit to `dev`: `refactor: Implement Problem/Physics abstraction layer`.
    * Merge `dev` to `master`: `milestone: Phase 2 complete. Core library abstracted.`.

---

## **Phase 3: Targeted Algorithmic Integration (DDM)**

**Objective:** Implement an advanced Domain Decomposition Method (DDM) solver, a key research frontier, within the new library structure.

### **Step 3.1: DDM Solver Skeleton and Test**

* **Planning:**
    * Implement an Overlapping Schwarz DDM preconditioner. This will be a new `SolverInterface` implementation.
    * This solver will manage local subdomain problems and MPI communication.
* **Test (`tests/test_solver_ddm_schwarz.cpp`):**
    * `test_ddm_schwarz_solver_solves_poisson`:
        1.  Run this test with MPI (`mpirun -np 4 ...`).
        2.  Use the Poisson MMS problem from Phase 1.
        3.  Instantiate `hpcfem::DdmSchwarzSolver`.
        4.  Instantiate `hpcfem::FemProblem` with this new solver.
        5.  `problem.solve();`
        6.  Assert that the *parallel* $L_2$ error norm is below `MMS_TOLERANCE`.
    * `test_ddm_schwarz_solver_converges_faster_than_jacobi`:
        1.  Run the same problem with a simple Jacobi solver (which is a basic DDM).
        2.  Run with the `DdmSchwarzSolver`.
        3.  Assert `ddm_iterations < jacobi_iterations`. (This proves the preconditioner is working).
* **Implementation:**
    * Create `include/hpcfem/solver_ddm_schwarz.hpp` and `src/solver_ddm_schwarz.cpp`.
    * `class DdmSchwarzSolver : public SolverInterface`.
    * The constructor will take the `mfem::ParMesh` and perform local setup (e.g., create local `HypreAmgSolver` for subdomain solves).
    * The `solve` method will implement the additive Schwarz iteration loop:
        1.  Update halo regions (MPI Send/Recv).
        2.  Solve local subdomain problems *in parallel*.
        3.  Check for global convergence.
    * Register new class in `docs/hpcfem-doc/`.
* **Test:** Run `ctest` with `mpirun`.
* **Refactor/Doc:**
    * Add Doxygen comments explaining the DDM algorithm and its MPI communication pattern.
    * Ensure no `TODO`s are left (Guideline #8).
    * Commit to `dev`: `feat: Implement overlapping Schwarz DDM solver`.

---

## **Phase 4: Rigorous Benchmarking & Analysis**

**Objective:** Quantitatively prove the new algorithm provides significant improvement over the baseline, as per document Phase 4.

### **Step 4.1: Benchmark Application**

* **Planning:**
    * Create a dedicated benchmark executable that can be controlled via command-line arguments.
* **Test (`tests/test_benchmark_cli.cpp`):**
    * `test_benchmark_cli_parses_solver`: Run `benchmark_app --solver amg` and check stdout for "Using AMG solver".
    * `test_benchmark_cli_parses_mesh`: Run `benchmark_app --mesh_order 2` and check stdout for "Using mesh order 2".
* **Implementation:**
    * Create `benchmarks/joule_heating_scaling/main.cpp`.
    * Implement command-line parsing (NO 3rd-party libs).
    * The `main` function will:
        1.  Parse args.
        2.  Load mesh.
        3.  Instantiate the correct `FemProblem`.
        4.  Select `hpcfem::HypreAmgSolver` or `hpcfem::DdmSchwarzSolver` based on args.
        5.  `problem.solve();`
        6.  Print timing, iteration count, and DOF count in a simple, parseable format (e.g., `DOFS: 100000\nTIME: 1.23\nITERS: 15\n`).
* **Test:** Run `ctest`.
* **Refactor/Doc:**
    * Add `benchmarks/joule_heating_scaling/README.md` explaining all CLI options.
    * Commit to `dev`: `feat: Create CLI benchmark application`.

### **Step 4.2: Scaling Study Automation**

* **Planning:**
    * Write a Python/Bash script to automate strong and weak scaling studies.
* **Test (Shell test in `tests/test_benchmarking_scripts.sh`):**
    * `test_run_strong_scaling_script`:
        1.  Run `benchmarks/run_scaling.py --type strong --min_procs 1 --max_procs 4`.
        2.  The script will execute `mpirun -np 1 ...`, `mpirun -np 2 ...`, etc.
        3.  It will parse the stdout and generate `results.csv`.
        4.  The test asserts that `results.csv` is created and contains 4 rows of data.
* **Implementation:**
    * Create `benchmarks/run_scaling.py` and `benchmarks/plot_results.py`.
    * The script will orchestrate the `mpirun` commands, collect output, and save it.
* **Test:** Run the shell test.
* **Refactor/Doc:**
    * Document the benchmarking workflow in `benchmarks/README.md`.
    * Commit to `dev`: `feat: Add benchmark automation scripts`.

### **Step 4.3: Stress Test and Analysis**

* **Planning:**
    * Execute the benchmark scripts for a large-scale problem (e.g., 1e7 DOFs) as the final stress test (Guideline #10).
* **Implementation (Execution):**
    * Run `benchmarks/run_scaling.py --type strong --procs 1,2,4,8,16,32 --dofs 1e7`.
    * Run `benchmarks/run_scaling.py --type weak --procs 1,2,4,8,16,32 --dofs_per_core 1e6`.
    * Run `benchmarks/plot_results.py` to generate `strong_scaling.png` and `weak_scaling.png`.
* **Test (Manual):**
    * Analyze the plots.
    * **Checklist:**
        * Does the DDM solver show better strong scaling than the monolithic AMG solver?
        * Does the DDM solver show good weak scaling (near-constant time)?
        * Did any run fail?
* **Refactor/Doc:**
    * If scaling is poor, stop and analyze (Guideline #11). The bottleneck is likely communication or poor subdomain solves. Refactor `DdmSchwarzSolver`.
    * Once scaling is confirmed, add the plots and a written analysis to `docs/benchmark_results_ddm.md`.
    * Commit to `dev`: `docs: Add DDM scaling results and analysis`.
    * Merge `dev` to `master`: `milestone: Phase 4 complete. DDM solver benchmarked and validated.`.

---

## **Phase 5: Advanced Optimization (Goal-Oriented AMR)**

**Objective:** Implement Goal-Oriented Adaptive Mesh Refinement (AMR) using the Dual Weighted Residual (DWR) method to demonstrate academic cutting-edge capability.

### **Step 5.1: Adjoint Problem Solver**

* **Planning:**
    * AMR via DWR requires solving an *adjoint* problem. For $L(u) = f$, the adjoint is $L^*(z) = g$, where $g$ is derived from the "quantity of interest" $J(u)$.
* **Test (`tests/test_amr_adjoint.cpp`):**
    * `test_adjoint_problem_assembly`:
        1.  Define a simple Poisson problem $L(u) = f$.
        2.  Define a goal $J(u) = \int_{\Omega} u \, dx$ (average value).
        3.  The test will instantiate `hpcfem::ElectrostaticsPhysics`.
        4.  Call a new method `physics->assembleAdjoint(J, A_adj, b_adj)`.
        5.  Assert that `A_adj` is the transpose of the original matrix $A$ and `b_adj` is the vector representation of $J$.
* **Implementation:**
    * Modify `physics_interface.hpp` to include `virtual void assembleAdjoint(...) = 0;`.
    * Implement this method in `physics_electrostatics.cpp`. This will involve using MFEM's `Transpose()` and `LinearForm` assembly for the goal $J$.
* **Test:** Run `ctest`.
* **Refactor/Doc:**
    * Document the new `assembleAdjoint` method and its mathematical basis.
    * Commit to `dev`: `feat: Add adjoint problem assembly to Physics layer`.

### **Step 5.2: DWR Error Estimator**

* **Planning:**
    * Implement the DWR error estimator, which computes local error indicators $\eta_K$ for each element $K$.
* **Test (`tests/test_amr_dwr_estimator.cpp`):**
    * `test_dwr_estimator_effectivity`:
        1.  Use the Poisson MMS problem on a coarse 4-element mesh.
        2.  Define a goal $J(u)$.
        3.  Solve the forward problem for $u_h$.
        4.  Solve the adjoint problem for $z_h$.
        5.  Create `hpcfem::DwrEstimator` and call `estimator.computeErrors(u_h, z_h)`.
        6.  Compute the true error $E = |J(u_{exact}) - J(u_h)|$.
        7.  Compute the total estimated error $\eta = \sum \eta_K$.
        8.  Assert that the effectivity index $\theta = \eta / E$ is close to 1.0 (e.g., $0.8 < \theta < 1.2$).
* **Implementation:**
    * Create `include/hpcfem/amr_dwr_estimator.hpp` and `src/amr_dwr_estimator.cpp`.
    * The `computeErrors` method will loop over all elements and compute the residual-based DWR formula.
    * Register new class in `docs/hpcfem-doc/`.
* **Test:** Run `ctest`.
* **Refactor/Doc:**
    * Document the DWR formula used.
    * Commit to `dev`: `feat: Implement DWR error estimator`.

### **Step 5.3: AMR Control Loop**

* **Planning:**
    * Create a top-level controller that runs the AMR loop: `SOLVE -> ESTIMATE -> MARK -> REFINE`.
* **Test (`tests/test_amr_controller.cpp`):**
    * `test_amr_controller_reduces_goal_error`:
        1.  Use the Poisson MMS problem with a known singularity (e.g., L-shaped domain).
        2.  Define a goal (e.g., average flux on one edge).
        3.  Create `hpcfem::AmrController`.
        4.  Run `controller.solveToTolerance(goal_tolerance)`.
        5.  **Checklist:**
            * Assert the loop terminates in $< \text{MAX_AMR_LOOPS}$.
            * Assert the final *goal error* is below `goal_tolerance`.
            * Assert the final mesh has *more* elements than the initial mesh.
            * Assert the refined elements are concentrated near the singularity.
* **Implementation:**
    * Create `include/hpcfem/amr_controller.hpp` and `src/amr_controller.cpp`.
    * The `solveToTolerance` method will contain the main loop.
    * Inside the loop:
        1.  `problem->solve();`
        2.  `adjoint_problem->solve();`
        3.  `estimator->computeErrors();`
        4.  `markElements();` (e.g., top 30% of elements).
        5.  `mesh->Refine(marked_elements);`
        6.  `problem->updateMesh(mesh);` (This requires re-assembling all physics).
* **Test:** Run `ctest`.
* **Refactor/Doc:**
    * This is a complex workflow. Refactor `FemProblem` to cleanly support mesh updates.
    * Add an example in `examples/amr_poisson/` to demonstrate the feature.
    * Commit to `dev`: `feat: Implement AMR controller with DWR`.
    * Merge `dev` to `master`: `milestone: Phase 5 complete. Goal-oriented AMR implemented.`.