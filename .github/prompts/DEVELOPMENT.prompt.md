---
mode: agent
---

# HPC-FEM Development Guidelines and Roadmap

Here is a comprehensive TDD roadmap to guide an AI developer in refactoring the existing project and sequentially building it into the state-of-the-art research framework described in your document. This roadmap adheres strictly to your provided "DEVELOPING GUIDELINES."

This roadmap is divided into five major phases, moving from foundational software architecture to cutting-edge algorithmic implementation.

## Phase 1: Foundation, Build System, and Core Abstractions

**Goal:** Establish a robust, modern, and extensible project layout that adheres to all development guidelines. This phase implements the vision from Part III of the provided document.

### Step 1.1: Modern CMake & Project Structure Refactor

* **PLAN:** Refactor the entire project to align with Part III and Guidelines 2, 14, 17, and 18. The current layout is a good start, but we must enforce strict target-based properties and dependency management.
* **TEST (Manual Check & Build):** This step is foundational. The "test" is to successfully build all targets (`hpcfem` library, all tests, all examples) on the `dev` branch after refactoring and to manually verify the new structure.
* **IMPLEMENT:**
    1.  **Top-level `CMakeLists.txt`:** Define the project, C++ standard (e.g., C++17), and enable testing. Add `hpcfem` as the core library target using `add_library(hpcfem)`.
    2.  **Dependency Management (Section 7.3):**
        * Integrate `CPM.cmake` (Section 7.2) for lightweight dependencies like GoogleTest.
        * Use `find_package` for heavy, system-level dependencies (MPI, MFEM, HYPRE). Assume these are provided by an environment manager like Spack (Section 7.3).
    3.  **`src/CMakeLists.txt`:**
        * Add sources to the `hpcfem` target.
        * Use `target_include_directories(hpcfem PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})` to make `src` the root include path. This enforces Guideline 18 (no `../`). All includes will be `<hpcfem/fem_problem.hpp>`.
        * Use `target_link_libraries(hpcfem PUBLIC mfem hypre)` to propagate dependencies.
    4.  **`tests/CMakeLists.txt`:**
        * Fetch GoogleTest using `CPMAddPackage`.
        * Define each test file as a separate executable (e.g., `add_executable(test_solver_interface test_solver_interface.cpp)`).
        * Link each test: `target_link_libraries(test_solver_interface PRIVATE hpcfem gtest_main)`.
        * Add tests to CTest: `add_test(NAME test_solver_interface COMMAND test_solver_interface)`.
    5.  **`example/` & `benchmark/`:** Refactor their `CMakeLists.txt` files identically to `tests/`.
    6.  **Code Cleanup:**
        * Move all `.hpp` files from `include/hpcfem` to `src/hpcfem` alongside their `.cpp` files (Guideline 17). Delete the `include` directory.
        * Update all `#include` paths in the project to reflect the new structure (e.g., `#include "fem_problem.hpp"` becomes `#include "hpcfem/fem_problem.hpp"`).
        * Place all code inside the `hpcfem` namespace (Guideline 2).
        * Refactor all names to match Guideline 16 (e.g., `FemProblem`, `solver_interface` -> `SolverInterface`).
        * Add Doxygen comments (Guideline 6) and update `docs/hpcfem-doc/naming_registry.md` (Guideline 5).

* **CHECKLIST:**
    * [ ] Project compiles fully with `cmake`, `cmake --build .`.
    * [ ] All tests run and pass with `ctest`.
    * [ ] `include` directory is deleted. `.hpp` and `.cpp` files are co-located in `src/hpcfem`.
    * [ ] All `#include` paths are relative to `src` (e.g., `#include "hpcfem/solver_interface.hpp"`).
    * [ ] All user code is in namespace `hpcfem`.
    * [ ] All naming conventions (Guideline 16) are applied.
    * [ ] No templates, lambdas, PCH, or macros are used (Guideline 1).

### Step 1.2: TDD - Refactor Physics and Solver Abstractions

* **PLAN:** The existing `PhysicsInterface` and `SolverInterface` are good, but they are not fully abstract or testable. We will solidify these interfaces.
* **TEST:**
    1.  **`test_solver_interface.cpp`:** Write a test `test_solver_interface_mock` that creates a mock `PhysicsInterface` and a mock `SolverInterface`. The test verifies that the `SolverInterface::solve` method is called once. This tests the *contract* of the interface.
    2.  **`test_physics_abstraction.cpp`:** Write a test `test_physics_assembly` that uses a mock `FemProblem` to verify that `PhysicsInterface::assemble` correctly calls the problem's assembly methods.
    3.  **`test_problem_abstraction.cpp`:** Update `test_problem_abstraction.cpp` to test the `FemProblem` class, verifying it can load a mesh and set up finite element spaces.
* **IMPLEMENT:**
    1.  `hpcfem/physics_interface.hpp`: Define `PhysicsInterface` with pure virtual methods: `virtual void assemble(FemProblem& problem) = 0;` and `virtual mfem::Operator& getOperator() = 0;`.
    2.  `hpcfem/solver_interface.hpp`: Define `SolverInterface` with a pure virtual method: `virtual void solve(mfem::Vector& x, mfem::Vector& b) = 0;` and `virtual void setOperator(mfem::Operator& op) = 0;`.
    3.  `hpcfem/fem_problem.hpp/cpp`: Refactor this class to *own* the mesh, finite element collection, and finite element space. It should provide public methods for physics classes to get references to these.
    4.  `hpcfem/physics_electrostatics.hpp/cpp`: Refactor this to be a concrete implementation of `PhysicsInterface`. It will *not* own the `FemProblem` but will take it as a reference in its `assemble` method.
    5.  `hpcfem/solver_hypre_amg.hpp/cpp`: Refactor this to be a concrete implementation of `SolverInterface`.
* **REFACTOR:** Ensure all classes from `fem_problem`, `physics_electrostatics`, and `solver_hypre_amg` now pass the new and existing tests.
* **DOCUMENT:** Update Doxygen comments for all refactored classes and register new names (Guideline 5).
* **CHECKLIST:**
    * [ ] All new and existing tests pass.
    * [ ] `PhysicsInterface` and `SolverInterface` contain only pure virtual methods (and a virtual destructor).
    * [ ] `PhysicsElectrostatics` and `SolverHypreAmg` are concrete, tested implementations.
    * [ ] `FemProblem` correctly manages mesh and FE space resources.

---

## Phase 2: Scalable Single-Physics Solvers (DDM + AMG)

**Goal:** Implement the scalable two-level Schwarz method with an AMG coarse solver, as described in Part I, Section 1. This directly addresses the performance bottleneck mentioned in the document.

### Step 2.1: TDD - One-Level Domain Decomposition (Benchmark for Failure)

* **PLAN:** Implement a basic, one-level overlapping Schwarz method as a preconditioner. The TDD goal here is to create a *scalability test* that *fails* (i.e., shows poor scaling), proving the document's point (Section 1.1).
* **TEST:**
    1.  **`tests/test_ddm_one_level.cpp`:** Write a test that solves a 3D Poisson problem (using `PhysicsElectrostatics`) on a 3D mesh (e.g., `testdata/testmesh_cube.mesh`) using MPI (e.g., 4 processes). The test passes if the solution converges and the result is correct.
    2.  **`benchmark/poisson_scaling_ddm1.cpp`:** Create a new benchmark that solves the 3D Poisson problem with increasing numbers of MPI processes (e.g., 2, 4, 8, 16). The benchmark *logs* the number of Krylov (GMRES/CG) iterations.
* **IMPLEMENT:**
    1.  `hpcfem/solver_one_level_schwarz.hpp/cpp`: Create a new class implementing `SolverInterface` (or rather, `mfem::Solver` to act as a preconditioner).
    2.  This class will use MFEM's parallel mesh partitioning (`mfem::ParMesh`) and finite element space (`mfem::ParFiniteElementSpace`).
    3.  The `solve` method will implement the additive Schwarz algorithm: parallel local solves on each subdomain's `A_i` (Section 1.2).
* **TESTING (Benchmark):** Run the `poisson_scaling_ddm1` benchmark. The "pass" condition for this step is observing and logging that the iteration count *increases significantly* as the process count grows. This confirms the "failure" of the one-level method (Section 1.1).
* **CHECKLIST:**
    * [ ] `test_ddm_one_level` passes (correctness).
    * [ ] `poisson_scaling_ddm1` benchmark runs and produces logs.
    * [ ] Logged data confirms that iteration count scales poorly with processor count.
    * [ ] New class `SolverOneLevelSchwarz` is documented and registered (Guideline 5).

### Step 2.2: TDD - Two-Level DDM with AMG Coarse Solver

* **PLAN:** Implement the two-level Schwarz method by adding a coarse-grid correction, as defined by the equation in Section 1.2. We will use our existing `SolverHypreAmg` as the coarse solver $A_H^{-1}$.
* **TEST:**
    1.  **`tests/test_ddm_two_level.cpp`:** Write a correctness test, similar to `test_ddm_one_level.cpp`, that solves the 3D Poisson problem using the new two-level preconditioner.
    2.  **`benchmark/poisson_scaling_ddm2.cpp`:** Copy the `poisson_scaling_ddm1` benchmark and modify it to use the new two-level solver.
* **IMPLEMENT:**
    1.  `hpcfem/solver_two_level_schwarz.hpp/cpp`: Create the new class.
    2.  It will require the same local solves as the one-level method (the $\sum R_i^T A_i^{-1} R_i$ term).
    3.  It will *also* construct the coarse grid. Use a standard Galerkin projection $A_H = \Psi^T A \Psi$ (Section 1.2). For the prolongation $\Psi$, use the built-in capabilities of `mfem::ParFiniteElementSpace` to define the coarse-to-fine mapping.
    4.  Instantiate `hpcfem::SolverHypreAmg` to solve the coarse problem $A_H^{-1}$ (Section 1.3).
    5.  Combine these in the preconditioner's `solve` method: $\tilde{B} = \Psi A_H^{-1} \Psi^T + \sum R_i^T A_i^{-1} R_i$ (Section 1.2).
* **TESTING (Benchmark):** Run the `poisson_scaling_ddm2` benchmark.
* **REFACTOR/VERIFY:** Compare the logged iteration counts from `poisson_scaling_ddm2` (two-level) against `poisson_scaling_ddm1` (one-level). The TDD pass condition is met if the iteration count for the two-level method remains *robustly low and nearly constant* as the number of processes increases. This demonstrates scalability (Section 1.2).
* **CHECKLIST:**
    * [ ] `test_ddm_two_level` passes (correctness).
    * [ ] `poisson_scaling_ddm2` benchmark runs and produces logs.
    * [ ] Logged data confirms that iteration count is low and (near) constant, proving scalability.
    * [ ] New class `SolverTwoLevelSchwarz` is documented and registered (Guideline 5).

---

## Phase 3: Physics-Informed Multiphysics Solvers

**Goal:** Implement block-structured preconditioners for coupled systems, moving from a "black-box" to a "physics-informed" approach as described in Part I, Section 2.

### Step 3.1: TDD - Define Coupled Electro-Thermal Problem

* **PLAN:** Implement the Joule Heating problem from Section 2.4. This requires a new "Thermal" physics class and a "Coupled" physics class.
* **TEST:**
    1.  **`tests/test_physics_thermal.cpp`:** Create a new test for a *standalone* transient heat problem (similar to `example/ex3_transient_heat.cpp`). This test verifies the thermal physics assembly.
    2.  **`tests/test_physics_coupling.cpp`:** Update this test (which is currently empty) to test the coupled electro-thermal problem. The test will:
        * Assemble the monolithic $2 \times 2$ block system (Section 2.4).
        * Solve it with a "black-box" direct solver (e.g., `mfem::UMFPackSolver`).
        * Verify the solution is correct (e.g., temperature increases where current flows).
* **IMPLEMENT:**
    1.  `hpcfem/physics_thermal.hpp/cpp`: Create a new `PhysicsThermal` class implementing `PhysicsInterface`.
    2.  `hpcfem/physics_joule_heating.hpp/cpp`: Create `PhysicsJouleHeating` implementing `PhysicsInterface`. This class will:
        * Internally hold instances of `PhysicsElectrostatics` and `PhysicsThermal`.
        * Override the `assemble` method to build the $2 \times 2$ block matrix $A = \begin{pmatrix} K_e(T) & 0 \\ C(V) & K_t \end{pmatrix}$ (Section 2.4) as an `mfem::BlockOperator`.
* **CHECKLIST:**
    * [ ] `test_physics_thermal` passes.
    * [ ] `test_physics_coupling` passes (i.e., the monolithic coupled problem can be assembled and solved correctly).
    * [ ] New classes `PhysicsThermal` and `PhysicsJouleHeating` are documented and registered (Guideline 5).

### Step 3.2: TDD - Block-Triangular (Gauss-Seidel) Preconditioner

* **PLAN:** Implement the block-triangular (Gauss-Seidel) preconditioner described in Sections 2.2 and 2.4. This is a physics-informed solver.
* **TEST:**
    1.  **`tests/test_solver_block_gauss_seidel.cpp`:** Create a new test. It will:
        * Set up the `PhysicsJouleHeating` problem.
        * Solve the system using GMRES preconditioned by the new `SolverBlockGaussSeidel`.
        * The test passes if the solver converges to the correct solution (comparing against the direct solver result from 3.1).
    2.  **`benchmark/joule_heating_solvers.cpp`:** Create a benchmark that solves the coupled problem with:
        * Monolithic AMG (e.g., `SolverHypreAmg` on the `mfem::BlockOperator`).
        * The new `SolverBlockGaussSeidel`.
        The benchmark logs iterations and solve time.
* **IMPLEMENT:**
    1.  `hpcfem/solver_block_gauss_seidel.hpp/cpp`: Create a new class implementing `mfem::Solver`.
    2.  The constructor will take references to the diagonal block solvers (e.g., two `SolverHypreAmg` instances, one for $K_e$ and one for $K_t$) (Section 2.4).
    3.  The `solve` method will implement the block forward substitution (Section 2.2):
        1.  Solve with $\tilde{K}_e$.
        2.  Update the right-hand side using $C(V)$.
        3.  Solve with $\tilde{K}_t$.
* **TESTING (Benchmark):** Run the `joule_heating_solvers` benchmark.
* **REFACTOR/VERIFY:** The TDD pass condition is met if the `SolverBlockGaussSeidel` converges in *significantly fewer iterations* and/or *less time* than the "black-box" monolithic AMG, proving the superiority of the physics-informed approach (Section 2.1).
* **CHECKLIST:**
    * [ ] `test_solver_block_gauss_seidel` passes (correctness).
    * [ ] `joule_heating_solvers` benchmark runs and produces logs.
    * [ ] Logs confirm the block solver is more efficient than the monolithic solver.
    * [ ] New class `SolverBlockGaussSeidel` is documented and registered (Guideline 5).

---

## Phase 4: Parallel-in-Time (PinT) Solvers

**Goal:** Implement a Parallel-in-Time (PinT) solver and demonstrate its effectiveness for parabolic problems and its failure/remedy for hyperbolic problems, as described in Part I, Section 3.

### Step 4.1: TDD - MGRIT for Parabolic Problems (Test for Success)

* **PLAN:** Implement the MGRIT (or Parareal) algorithm (Section 3.1) and apply it to the parabolic (diffusive) transient heat problem.
* **TEST:**
    1.  **`tests/test_solver_mgrit_parabolic.cpp`:**
        * Set up the `PhysicsTransientHeat` problem.
        * First, solve it using standard, sequential time-stepping (e.g., Backward Euler) and save the final solution vector.
        * Second, solve the *same* problem using the new MGRIT solver.
        * The test passes if the final solution from MGRIT *exactly matches* the sequential solution after a few iterations.
* **IMPLEMENT:**
    1.  `hpcfem/solver_mgrit.hpp/cpp`: Create a new class. This is complex. It will need to:
        * Implement the MGRIT V-cycle (Section 3.1).
        * Define a "fine" propagator $\mathcal{F}$ (e.g., one step of Backward Euler).
        * Define a "coarse" propagator $\mathcal{G}$ (e.g., one step of Backward Euler with a much larger time step $\Delta T = N \Delta t$) (Section 3.1).
        * The "all-at-once" system is solved iteratively. This requires MPI parallelism across time-steps.
* **CHECKLIST:**
    * [ ] `test_solver_mgrit_parabolic` passes (correctness).
    * [ ] The MGRIT solver converges in a few iterations (e.S., 3-5) for the parabolic problem.
    * [ ] New class `SolverMGRIT` is documented and registered (Guideline 5).

### Step 4.2: TDD - MGRIT for Hyperbolic Problems (Test for Failure & Fix)

* **PLAN:** Demonstrate the failure of standard PinT for hyperbolic problems (Section 3.2) and then implement the physics-aware fix (Section 3.3).
* **TEST (Cycle):**
    1.  **`tests/test_physics_advection.cpp`:** Create a new test for a simple 3D advection equation. Test passes if the sequential time-stepped solution is correct.
    2.  **`tests/test_solver_mgrit_hyperbolic_fail.cpp`:** Run the *exact same* `SolverMGRIT` from Step 4.1 on the new `PhysicsAdvection` problem. The test *passes* if the solver *fails to converge* in a reasonable number of iterations (e.g., > 50), confirming the "hyperbolic challenge" (Section 3.2).
    3.  **`tests/test_solver_mgrit_hyperbolic_pass.cpp`:** This is the final test. It runs a *modified* MGRIT solver on the advection problem. The test passes if it *now converges rapidly* (e.g., < 10 iterations).
* **IMPLEMENT:**
    1.  `hpcfem/physics_advection.hpp/cpp`: Create a new `PhysicsAdvection` class.
    2.  Run the test `test_solver_mgrit_hyperbolic_fail` and verify it fails to converge.
    3.  **Refactor `SolverMGRIT`:**
        * Modify the `SolverMGRIT` class to accept a *strategy* or *sub-class* for the coarse-grid operator $\mathcal{G}$.
        * Create `CoarseOpStandard` (the original, failing $\Delta T$ version).
        * Create `CoarseOpSemiLagrangian` (the fix). This operator traces the characteristic curves backward in time to build the coarse propagator (Section 3.3).
* **TESTING:** Run `test_solver_mgrit_hyperbolic_pass` (using `CoarseOpSemiLagrangian`). Verify it passes.
* **CHECKLIST:**
    * [ ] `test_physics_advection` passes.
    * [ ] `test_solver_mgrit_hyperbolic_fail` *logs* a failure to converge.
    * [ ] `test_solver_mgrit_hyperbolic_pass` *passes* with the new characteristic-based operator.
    * [ ] `SolverMGRIT` is refactored, and new classes (`PhysicsAdvection`, `CoarseOpSemiLagrangian`) are documented and registered (Guideline 5).

---

## Phase 5: AI-Enhanced Numerical Solvers

**Goal:** Integrate machine learning to learn and optimize solver components, as described in Part II. This requires a C++/Python bridge.

### Step 5.1: TDD - GNN-Based AMG Prolongation Operator

* **PLAN:** Use a Graph Neural Network (GNN) to learn the AMG prolongation operator $P$, replacing the standard heuristic in `SolverHypreAmg` (or a new, custom AMG) (Section 4).
* **TEST:**
    1.  **`tests/test_python_bridge.cpp`:** Create a test to verify that C++ can initialize the Python interpreter (using `hpc-fem-playground` conda env, Guideline 19), import a simple Python module (e.g., `torch`), and call a function.
    2.  **`tests/test_gnn_amg.cpp`:** This is the key test.
        * Set up a "hard" problem (e.g., `PhysicsElectrostatics` with highly anisotropic coefficients).
        * Solve it using the standard `SolverHypreAmg` and log the iteration count (this is our baseline).
        * Solve it using the new `SolverGnnAmg`.
        * The test passes if the iteration count for `SolverGnnAmg` is *less than or equal to* the baseline.
* **IMPLEMENT:**
    1.  **Python Bridge:** Integrate a C++/Python interop library (e.g., `pybind11`). Configure CMake to build it.
    2.  **Python GNN:** In `scripts/` (or a new `ml/` dir), create a Python script `gnn_amg.py` that uses `pytorch_geometric` to define a GNN (Section 4.2).
    3.  **Python Training:** Create a training script (`train_gnn.py`) that loads matrices, trains the GNN using the *unsupervised energy-minimizing loss function* $\mathcal{L}(P) = \sum_{k} \| E \cdot e_{smooth}^{(k)} \|_A^2$ (Section 4.3), and saves the trained model.
    4.  `hpcfem/solver_gnn_amg.hpp/cpp`: Create a new `SolverInterface`. This class will:
        * Call the Python bridge to load the pre-trained GNN model.
        * During setup, extract the matrix graph, pass it to the GNN, and get back the prolongation operator $P$.
        * Build the rest of the AMG hierarchy using this $P$.
* **CHECKLIST:**
    * [ ] `test_python_bridge` passes.
    * [ ] `train_gnn.py` script runs and saves a model file.
    * [ ] `test_gnn_amg` passes (i.e., GNN-based solver is at least as good as the baseline).
    * [ ] New class `SolverGnnAmg` is documented and registered (Guideline 5).

### Step 5.2: TDD - RL-Based Hyperparameter Tuning

* **PLAN:** Use Reinforcement Learning (RL) to *dynamically* tune a solver parameter online, as described in Section 5. We will target the relaxation factor $\omega$ in a standard iterative solver (e.g., SOR) or the penalty parameter $\rho$ in ADMM (Section 5.2).
* **TEST:**
    1.  **`benchmark/rl_tuning_vs_static.cpp`:** Create a benchmark to solve a specific problem (e.g., `PhysicsElectrostatics`).
        * **Baseline:** Run the solver (e.g., `mfem::SOR`) with the *theoretically optimal static* relaxation parameter $\omega$. Log the wall-clock time.
        * **RL-Tuned:** Run the new `SolverRlTunedSor`.
        * The TDD pass condition is that the *total wall-clock time* for the RL-tuned solver (including RL agent overhead) is *less than* the baseline static solver's time.
* **IMPLEMENT:**
    1.  **Python RL Agent:** In `ml/`, create `rl_agent.py`. This defines the RL agent (e.g., a Deep Q-Network) that implements the policy $\pi(a|s)$ (Section 5.3).
    2.  **Training:** Create `train_rl.py` to train the agent by repeatedly running the solver, feeding it `State` (residual history) and `Reward` (residual reduction) (Section 5.2).
    3.  `hpcfem/solver_rl_tuned_sor.hpp/cpp`: Create a new `SolverInterface`. This C++ class will:
        * Call the Python bridge to load the pre-trained RL agent.
        * In its `solve` loop, at each iteration $k$:
            1.  Get the current `State` $s_k$ (residual norm, etc.).
            2.  Call the RL agent $\pi(a_k|s_k)$ to get the *action* $a_k$ (the new parameter $\omega_k$).
            3.  Perform the SOR iteration using $\omega_k$.
            4.  Calculate the `Reward` $r_k$ and (if training) send it back to the agent.
* **CHECKLIST:**
    * [ ] `train_rl.py` script runs and saves a trained agent model.
    * [ ] `benchmark/rl_tuning_vs_static.cpp` benchmark runs.
    * [ ] Benchmark logs confirm the RL-tuned solver is *faster* in wall time than the optimal static solver, proving the value of online adaptation (Section 5.3).
    * [ ] New class `SolverRlTunedSor` is documented and registered (Guideline 5).