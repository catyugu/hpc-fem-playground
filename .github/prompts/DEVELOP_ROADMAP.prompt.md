---
mode: agent
---

Here is a comprehensive, multi-phase roadmap designed to guide an AI coder from a baseline implementation to a cutting-edge, academic-level optimization, as defined by your research guide.

This roadmap follows the phased development structure outlined in your document and integrates the advanced algorithmic frontiers from Part II as the ultimate research goals. Each phase includes strict implementation hints and automated checklists.

---

## **Phase 1: Baseline & Replication**

**Objective:** Establish a correct, validated, and parallel starting point for a representative multiphysics problem. This baseline is the benchmark against which all future optimizations will be measured.

### **Implementation Hints**

1.  **Problem Selection:** Start with a canonical coupled problem, such as **Joule heating** (electro-thermal). This involves coupling a current/potential problem (like MFEM's `ex1`) with a thermal problem (also `ex1`).
2.  **Benchmark Identification:** Find a well-regarded academic paper or benchmark study on this problem that includes published results (e.g., a convergence plot, max temperature for a given setup).
3.  **Parallel-First:** Implement the baseline *from the beginning* using MFEM's parallel classes: `ParMesh`, `ParFiniteElementSpace`, `ParGridFunction`, `HypreParMatrix`. Do not start with serial code.
4.  **Baseline Solver:** Use the standard, robust solver setup: a preconditioned conjugate gradient (PCG) solver preconditioned with **HYPRE's BoomerAMG**. This is your `Solver_v1`.
5.  **Discretization:** Use MFEM's high-order element capabilities. Start with $H^1$ elements of order $p=2$ or $p=3$. This aligns with MFEM's core design philosophy of increasing arithmetic intensity.
6.  **I/O:** Ensure the simulation saves its final state (mesh and solution `ParGridFunction`) and logs key metrics (e.g., "Number of PCG iterations", "Total solver time", "Max temperature").

### **Automation Checklist**

```[ ]``` **Task:** Identify and select a suitable MFEM miniapp (e.g., `ex1p`, `ex11p`) as a starting point.
```[ ]``` **Task:** Load a 3D mesh relevant to an IC feature (e.g., a TSV or interconnect).
```[ ]``` **Task:** Implement the coupled Joule heating problem by defining the necessary `ParBilinearForm` and `ParLinearForm` objects.
```[ ]``` **Task:** Configure the solver to use `HypreBoomerAMG` as a preconditioner for `PCG`.
```[ ]``` **Task:** Run the parallel simulation (e.g., `mpirun -np 4 ./your_app`).
```[ ]``` **Task:** Post-process the output and extract the primary validation metric.
```[ ]``` **Validation:** **The primary metric (e.g., "Max Temperature") matches the published benchmark result within a 1% tolerance**.
```[ ]``` **Version Control:** All code from this phase is committed to a `git` repository with a tag `v1.0-baseline`.

---

## **Phase 2: Encapsulation & Abstraction**

**Objective:** Transform the Phase 1 script into a modular, robust, and extensible software library. This is a critical software engineering step that enables rapid and reliable research.

### **Implementation Hints**

1.  **Refactor (Do Not Modify MFEM):** Adhere to the cardinal rule: **Use, Don't Modify**. All new code will be *on top* of MFEM, using its public API.
2.  **Solver Abstraction:** Create an abstract base class `SolverInterface` with a pure virtual method: `virtual void Solve(const Vector& b, Vector& x) = 0;`. This decouples the problem from the solver.
3.  **Concrete Baseline Solver:** Create a `BaselineAMGSolver` class that inherits from `SolverInterface`. Move the PCG+BoomerAMG logic from Phase 1 into this class's `Solve()` method.
4.  **Physics Abstraction:** Create a `PhysicsProblem` class. This class will be responsible for owning the `ParMesh` and `ParFiniteElementSpace` and will have methods like `AssembleSystem(BilinearForm& a, LinearForm& b)`.
5.  **Application Layer:** The `main.cpp` (Application Layer) should now be very clean. It should instantiate the `PhysicsProblem` and the `SolverInterface` (using `std::unique_ptr`) and be able to switch solvers with a single line of code.
6.  **Testing:** Implement unit tests (e.g., using Google Test) for any new utility functions (e.g., a custom boundary condition class).

### **Automation Checklist**

```[ ]``` **Task:** Create `SolverInterface.hpp`.
```[ ]``` **Task:** Create `BaselineAMGSolver.hpp/.cpp` inheriting from `SolverInterface` and encapsulating Phase 1 logic.
```[ ]``` **Task:** Create `PhysicsProblem.hpp/.cpp` encapsulating the problem setup (mesh, FE space, forms).
```[ ]``` **Task:** Refactor `main.cpp` to use these new classes and C++ best practices (e.g., RAII, `std::unique_ptr`).
```[ ]``` **Task:** Add Doxygen-style API documentation to all new public classes and methods.
```[ ]``` **Task:** Create a simple unit test for one component of the new library.
```[ ]``` **Validation:** **The refactored code compiles and, when run, produces bit-for-bit identical numerical results to the Phase 1 baseline.**
```[ ]``` **Version Control:** All refactored code is committed with a tag `v2.0-abstracted`.

---

## **Phase 3: Targeted Algorithmic Integration (The Research Frontier)**

**Objective:** Implement one of the advanced, cutting-edge algorithms from Part II within your new library structure. This new algorithm will inherit from `SolverInterface` for a direct, fair comparison with the baseline.

**Select one of the following paths (3A, 3B, or 3C) based on your research goals.**

### **Path 3A: Domain Decomposition Method (DDM)**

**Goal:** Implement a scalable DDM, such as an overlapping Schwarz preconditioner, which is highly suited for multiphysics and parallel scaling.

* **Implementation Hints:**
    1.  Create a new class `DDMSchwarzSolver` that inherits from `SolverInterface`.
    2.  This is a "divide and conquer" method. Your `ParMesh` already provides the domain partitioning.
    3.  The core iterative loop involves:
        * **Local Solves:** Solving the problem on each subdomain. This can re-use your `BaselineAMGSolver` for the local problems!
        * **Halo Exchange:** Exchanging solution data at the subdomain interfaces (halos) using MPI. Leverage `ParGridFunction::ExchangeFaceNbrData()`.
    4.  A key research area is the *transmission condition* (the boundary condition applied at the interfaces). Start with simple Dirichlet conditions.

* **Automation Checklist (3A):**
    ```[ ]``` **Task:** Create `DDMSchwarzSolver.hpp/.cpp`.
    ```[ ]``` **Task:** Implement the local subdomain solver logic.
    ```[ ]``` **Task:** Implement the interface data exchange (halo exchange) logic.
    ```[ ]``` **Task:** Implement the iterative Schwarz procedure as the `Solve()` method.
    ```[ ]``` **Task:** Add a unit test for the DDM solver on a simple 2D Poisson problem.
    ```[ ]``` **Validation:** **The `DDMSchwarzSolver` converges to the same solution as the `BaselineAMGSolver` (within solver tolerance).**
    ```[ ]``` **Version Control:** All DDM code is committed with a tag `v3.0-ddm`.

### **Path 3B: Goal-Oriented Adaptive Mesh Refinement (AMR)**

**Goal:** Implement a goal-oriented $h$-AMR strategy based on the Dual Weighted Residual (DWR) method to optimize the mesh for a specific engineering quantity.

* **Implementation Hints:**
    1.  **Goal:** Define a "quantity of interest," e.g., the maximum temperature at a specific hotspot location.
    2.  **Dual Problem:** Formulate the *adjoint* (dual) PDE associated with your quantity of interest.
    3.  **Implement Adjoint Solver:** Create a new solver (e.g., `AdjointSolver`) to solve this dual problem.
    4.  **Error Indicator:** Implement the DWR error estimator. This involves integrating the *primal residual* against the *dual solution* element-by-element.
    5.  **Refinement Loop:** Use MFEM's built-in AMR capabilities (e.g., `ParMesh::Refine()`). Feed your DWR indicators into an `ThresholdRefiner` to guide $h$-refinement.

* **Automation Checklist (3B):**
    ```[ ]``` **Task:** Formulate the adjoint PDE for your quantity of interest.
    ```[ ]``` **Task:** Implement the `AdjointSolver` and the `BilinearForm` for the adjoint problem.
    ```[ ]``` **Task:** Implement the DWR error indicator computation.
    ```[ ]``` **Task:** Create the main AMR loop (Solve Primal -> Solve Dual -> Compute Indicators -> Refine Mesh).
    ```[ ]``` **Task:** Add a regression test for the AMR solver.
    ```[ ]``` **Validation:** **Demonstrate that the AMR-generated mesh achieves the same accuracy for the *quantity of interest* using significantly fewer (e.g., >50% less) degrees of freedom than a uniformly refined mesh.**
    ```[ ]``` **Version Control:** All AMR code is committed with a tag `v3.0-amr`.

### **Path 3C: Model Order Reduction (MOR)**

**Goal:** Implement a data-driven, projection-based MOR (e.g., POD) to enable near-real-time solutions for "many-query" tasks like parameter sweeps.

* **Implementation Hints:**
    1.  **Offline Stage:**
        * Run your high-fidelity (Phase 2) model for several different parameter values (e.g., 20 different material conductivities).
        * Store all solution vectors (`ParGridFunction`) as "snapshots".
        * Perform Singular Value Decomposition (SVD) on the snapshot matrix to extract the $r$ most dominant POD modes (the reduced basis).
    2.  **Online Stage:**
        * Create a `ReducedOrderSolver` class.
        * For a *new* parameter, project the full-system operators (`BilinearForm`, `LinearForm`) onto the low-dimensional POD basis. This creates a tiny $r \times r$ dense system.
        * Solve this small system and project the solution back to the full-dimensional space.

* **Automation Checklist (3C):**
    ```[ ]``` **Task:** Implement an "offline" script that runs the Phase 2 model in a loop and saves snapshots.
    ```[ ]``` **Task:** Implement a script (e.g., in Python/SciPy or using a C++ linear algebra library) to perform SVD and generate the POD basis.
    ```[ ]``` **Task:** Create `ReducedOrderSolver.hpp/.cpp`.
    ```[ ]``` **Task:** Implement the operator projection logic in the `ReducedOrderSolver`.
    ```[ ]``` **Validation (Accuracy):** **For a new test parameter, the MOR solution has a <1% error compared to the full Phase 2 simulation.**
    ```[ ]``` **Validation (Speed):** **The "online" solve time is orders of magnitude faster (e.g., >1000x) than the full Phase 2 simulation.**
    ```[ ]``` **Version Control:** All MOR code is committed with a tag `v3.0-mor`.

---

## **Phase 4: Rigorous Benchmarking & Analysis**

**Objective:** Quantitatively prove that your Phase 3 algorithm provides a significant, measurable improvement over the Phase 2 baseline, and generate the final results for your research.

### **Implementation Hints**

1.  **Automation:** This phase is all about automation. Create shell scripts (`.sh`) or Python scripts to run your experiments. These scripts must:
    * Compile the code.
    * Execute the simulation with varying parameters (e.g., core count, problem size, refinement level).
    * Parse the log files (e.g., using `grep` or `awk`) to extract key metrics (Time, Iterations, DOFs, Error).
2.  **Strong Scaling Study:**
    * **Goal:** Measure speedup for a *fixed total problem size*.
    * **Action:** Run your Phase 3 solver on $N=1, 2, 4, 8, 16, 32...$ cores.
    * **Plot:** Runtime vs. Number of Cores.
3.  **Weak Scaling Study:**
    * **Goal:** Measure performance as the problem and machine size grow together (fixed *problem size per core*).
    * **Action:** Run on $N=1$ core (problem size $S$), $N=2$ cores (problem size $2S$), $N=4$ cores (problem size $4S$), etc.
    * **Plot:** Runtime vs. Number of Cores (Ideal runtime is a flat horizontal line).
4.  **Data Analysis:** Use Python (Pandas, Matplotlib) to process the raw output data and generate the final, publication-quality plots.
5.  **Profiling:** Use a profiler (e.g., `gprof`, `valgrind`, or NVIDIA's `nsight`) to identify the *next* bottleneck in your new code.

### **Automation Checklist**

```[ ]``` **Task:** Create an automation script (e.g., `run_strong_scaling.sh`) that runs the Phase 2 baseline and Phase 3 algorithm from 1 to 64 cores.
```[ ]``` **Task:** Create an automation script (e.g., `run_weak_scaling.sh`) for both solvers.
```[ ]``` **Task:** Create a Python script (e.g., `plot_results.py`) to read all log files and automatically generate scaling plots comparing the two solvers.
```[ ]``` **Task:** Run a profiler on the Phase 3 solver at high core count and identify the top 3 most expensive functions.
```[ ]``` **Validation (Final Report):** **Produce a plot and a data table that clearly demonstrates a quantifiable, significant improvement** (e.g., "The Phase 3 DDM solver achieves 5x speedup over the baseline AMG solver for a 1M DOF problem" or "The Phase 3 AMR solver achieves the target accuracy with 80% fewer DOFs").