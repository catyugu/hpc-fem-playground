---
mode: agent
---

# HPC-FEM Development Guidelines and Roadmap

Here is a comprehensive TDD roadmap to refactor the `hpc-fem-playground` to support user-defined FE spaces and integrators, strictly adhering to your developing guidelines.

**New Git Branch:** `feature/custom_fem_interfaces` (Guideline 9)
**Primary Documentation File:** `docs/hpcfem-doc/custom_fem_interfaces.md` (Create and update this file throughout all steps) (Guideline 22)

-----

### Part 1: Refactor Physics Modules for FE-Space Injection

**Goal:** Refactor `PhysicsInterface` and its concrete implementations (e.g., `PhysicsElectrostatics`, `PhysicsThermal`) to *receive* a `mfem::ParFiniteElementSpace` rather than creating one internally. This allows the user to construct and provide any `mfem::ParFiniteElementSpace` (including custom derived classes) at the application level.

#### Step 1.1: Planning

  * **Analysis:** `PhysicsElectrostatics` and `PhysicsThermal` currently create their own `fec` and `fespace` in their constructors. The `setupOn` method then redundantly checks this internal `fespace` against the one provided by `FemProblem`.
  * **Refactor Plan:**
    1.  Modify `PhysicsInterface` to add a protected member: `mfem::ParFiniteElementSpace* fespace = nullptr;`.
    2.  Remove the `fec` and `fespace` members from `PhysicsElectrostatics` and `PhysicsThermal`.
    3.  Change the constructors of `PhysicsElectrostatics` and `PhysicsThermal` to no longer accept an `int order`. They only need the `mfem::ParMesh& pmesh` to query boundary attributes.
    4.  Modify the `setupOn(mfem::ParFiniteElementSpace &fespace_in)` method in all `Physics` children:
          * Remove the redundant `fespace->GetVDim() != fespace_in.GetVDim()` check.
          * Assign the `fespace` from the `PhysicsInterface` base class: `this->fespace = &fespace_in;`.
          * All internal forms (e.g., `a`, `b`) must now be created using `this->fespace`.
  * **Affected Files:**
      * `src/hpcfem/physics_interface.hpp`
      * `src/hpcfem/physics_electrostatics.hpp` / `.cpp`
      * `src/hpcfem/physics_thermal.hpp` / `.cpp`
      * All files in `example/` and `tests/` that instantiate these physics modules.

#### Step 1.2: Write Failing Tests

1.  **Modify `tests/test_problem_abstraction.cpp`:**
      * **Hint:** This test will fail to compile. Update the instantiation of `PhysicsElectrostatics` to use the new constructor: `auto physics = std::make_unique<hpcfem::PhysicsElectrostatics>(*pmesh);`.
      * This test now fails until the implementation in Step 1.3 is complete.
2.  **Create New Test: `tests/test_custom_fespace_injection.cpp`**
      * **Checklist:**
          * [ ] Initialize MPI and load `testdata/testmesh_cube.mesh` as a `mfem::ParMesh`.
          * [ ] Create a *non-default* FE collection, e.g., `auto *fec = new mfem::L2_FECollection(order, pmesh.Dimension());`.
          * [ ] Create the `mfem::ParFiniteElementSpace* fespace` using this `fec`.
          * [ ] Create `auto problem = std::make_unique<hpcfem::FemProblem>(pmesh, fespace);`.
          * [ ] Create `auto physics = std::make_unique<hpcfem::PhysicsElectrostatics>(*pmesh);`.
          * [ ] Call `problem->addPhysics(physics.get());`.
          * [ ] **Assertion:** The test must `ASSERT(physics->getBilinearForm() != nullptr)` (after `setup()` is called by `addPhysics`).
          * [ ] **Assertion:** The test must `ASSERT(physics->getBilinearForm()->FESpace() == fespace)`. This proves the `physics` module is using the externally-created L2 space, not one it created internally.
      * This test will fail to compile or run until Step 1.3 is complete.

#### Step 1.3: Implementation

  * **Checklist:**
      * [ ] **`src/hpcfem/physics_interface.hpp`:** Add `protected: mfem::ParFiniteElementSpace* fespace = nullptr;`.
      * [ ] **`src/hpcfem/physics_electrostatics.hpp`:**
          * Change constructor to `explicit PhysicsElectrostatics(mfem::ParMesh &pmesh);`.
          * Remove `fec` and `fespace` members.
      * [ ] **`src/hpcfem/physics_electrostatics.cpp`:**
          * Update constructor to remove `fec` and `fespace` creation.
          * Update `setupOn` as described in "Planning" (remove check, set `this->fespace`, use `this->fespace` for `a` and `b`).
      * [ ] **`src/hpcfem/physics_thermal.hpp / .cpp`:** Apply the exact same refactoring as for `PhysicsElectrostatics`.

#### Step 1.4: Testing, Refactoring, and Documentation

1.  **Testing:**
      * **Checklist:**
          * [ ] Run `tests/test_problem_abstraction.cpp`. It must pass.
          * [ ] Run `tests/test_custom_fespace_injection.cpp`. It must pass.
2.  **Refactoring (Guideline 7):**
      * **Checklist:**
          * [ ] Update `example/ex0_cylinder_dirichlet.cpp`.
          * [ ] Update `example/ex1_cube_mixed_bc.cpp`.
          * [ ] Update `example/ex3_transient_heat.cpp`.
          * **Hint:** The only change in these files should be removing the `order` argument from the `Physics...` constructor (e.g., `auto physics = std::make_unique<hpcfem::PhysicsElectrostatics>(*pmesh, order);` becomes `auto physics = std::make_unique<hpcfem::PhysicsElectrostatics>(*pmesh);`).
          * [ ] Verify all examples compile and run correctly.
3.  **Documentation (Guideline 6, 22):**
      * **Checklist:**
          * [ ] Update Doxygen comments for all modified constructors and methods.
          * [ ] Update `docs/hpcfem-doc/custom_fem_interfaces.md` to mark Part 1 as complete.

-----

### Part 2: Refactor Physics Modules for Integrator Injection

**Goal:** Refactor `PhysicsInterface` and children to *accept* user-created integrators, rather than hard-coding them in `setupOn`. This allows the user to provide standard MFEM integrators or their own custom-derived classes.

#### Step 2.1: Planning

  * **Analysis:** `Physics...::setupOn` methods hard-code the creation of coefficients (e.g., `epsilonCoeff`) and integrators (e.g., `mfem::DiffusionIntegrator`). This logic must be removed from the `Physics` modules and moved to the user-facing application (`example/` or `test/`).
  * **Refactor Plan:**
    1.  Create a new header `src/hpcfem/integrator_info.hpp` (Guideline 14, 24). This file will define structs to hold integrator pointers and their optional markers. (Guideline 4: No nested classes).
    2.  **`src/hpcfem/integrator_info.hpp` content:**
        ```cpp
        #include "mfem.hpp"
        #include <memory>

        namespace hpcfem {

        // Info struct for Bilinear Domain Integrators
        struct BilinearIntegratorInfo {
            std::unique_ptr<mfem::BilinearFormIntegrator> integrator;
            std::unique_ptr<mfem::Array<int>> marker; // nullptr for all domains
        };

        // Info struct for Linear Domain Integrators
        struct LinearIntegratorInfo {
            std::unique_ptr<mfem::LinearFormIntegrator> integrator;
            std::unique_ptr<mfem::Array<int>> marker; // nullptr for all domains
        };

        // Info struct for Bilinear Boundary Integrators
        struct BilinearBoundaryIntegratorInfo {
            std::unique_ptr<mfem::BilinearFormIntegrator> integrator;
            std::unique_ptr<mfem::Array<int>> marker; // MUST NOT be null
        };

        // Info struct for Linear Boundary Integrators
        struct LinearBoundaryIntegratorInfo {
            std::unique_ptr<mfem::LinearFormIntegrator> integrator;
            std::unique_ptr<mfem::Array<int>> marker; // MUST NOT be null
        };

        } // namespace hpcfem
        ```
    3.  Modify `src/hpcfem/physics_interface.hpp`:
          * `#include "hpcfem/integrator_info.hpp"`.
          * Add `protected` vectors to store the integrators:
              * `std::vector<hpcfem::BilinearIntegratorInfo> bilinearIntegrators;`
              * `std::vector<hpcfem::LinearIntegratorInfo> linearIntegrators;`
              * `std::vector<hpcfem::BilinearBoundaryIntegratorInfo> bilinearBoundaryIntegrators;`
              * `std::vector<hpcfem::LinearBoundaryIntegratorInfo> linearBoundaryIntegrators;`
          * Add `public virtual` methods to add integrators:
              * `virtual void addBilinearIntegrator(std::unique_ptr<mfem::BilinearFormIntegrator> integ, std::unique_ptr<mfem::Array<int>> marker = nullptr);`
              * `virtual void addLinearIntegrator(std::unique_ptr<mfem::LinearFormIntegrator> integ, std::unique_ptr<mfem::Array<int>> marker = nullptr);`
              * `virtual void addBilinearBoundaryIntegrator(std::unique_ptr<mfem::BilinearFormIntegrator> integ, std::unique_ptr<mfem::Array<int>> marker);`
              * `virtual void addLinearBoundaryIntegrator(std::unique_ptr<mfem::LinearFormIntegrator> integ, std::unique_ptr<mfem::Array<int>> marker);`
    4.  Modify `src/hpcfem/physics_interface.cpp`: Implement these `add...` methods. They just construct the appropriate `...Info` struct and `std::move` the `unique_ptr`s into the correct `std::vector`.
    5.  Modify `PhysicsElectrostatics` (and `PhysicsThermal`):
          * **`hpp`:** Remove all `set...` methods (`setEpsilon`, `setSource`, `setNeumann`, `setDirichlet`). Remove all `mfem::Coefficient*`, `mfem::Vector*`, and marker `mfem::Array<int>` members.
          * **`cpp`:**
              * Update constructor: Remove all coefficient and marker initialization.
              * Update `setupOn`: Remove *all* `new ...Integrator` and `new ...Coefficient` lines.
              * **Implement New Logic in `setupOn`:**
                ```cpp
                // This logic MUST follow the Part 1 refactor
                delete a; a = new mfem::ParBilinearForm(this->fespace);
                delete b; b = new mfem::ParLinearForm(this->fespace);

                // Add domain bilinear integrators
                for (auto &info : bilinearIntegrators) {
                    if (info.marker == nullptr) {
                        a->AddDomainIntegrator(info.integrator.release()); // MFEM takes ownership
                    } else {
                        a->AddDomainIntegrator(info.integrator.release(), *info.marker);
                    }
                }
                bilinearIntegrators.clear(); // Clear vector of now-null ptrs

                // ... Repeat for linearIntegrators on b ...
                linearIntegrators.clear();

                // ... Repeat for bilinearBoundaryIntegrators on a ...
                bilinearBoundaryIntegrators.clear();

                // ... Repeat for linearBoundaryIntegrators on b ...
                linearBoundaryIntegrators.clear();
                ```
                **Note:** This implementation (`release()` and `clear()`) means `setupOn` is a one-time operation. This is correct and MFEM-compliant, as the forms take ownership of the integrator pointers.

#### Step 2.2: Write Failing Tests

1.  **Modify `tests/test_problem_abstraction.cpp`:**
      * **Hint:** This test will fail to compile as all `set...` methods are gone.
      * **New Test Logic:**
          * [ ] `auto physics = std::make_unique<hpcfem::PhysicsElectrostatics>(*pmesh);`
          * [ ] Create coefficients: `auto *epsCoeff = new mfem::ConstantCoefficient(1.0);`
          * [ ] Create integrators: `auto diffInteg = std::make_unique<mfem::DiffusionIntegrator>(*epsCoeff);`
          * [ ] Add integrators: `physics->addBilinearIntegrator(std::move(diffInteg));`
          * [ ] **(Repeat for a source term)** `auto *srcCoeff = new mfem::ConstantCoefficient(1.0);`
          * [ ] `auto srcInteg = std::make_unique<mfem::DomainLFIntegrator>(*srcCoeff);`
          * [ ] `physics->addLinearIntegrator(std::move(srcInteg));`
          * [ ] Add `physics` to `problem`, `setup`, `solve`, and assert a non-zero solution.
      * This test fails until Step 2.3 is complete.
2.  **Create New Test: `tests/test_custom_integrator.cpp`**
      * **Hint:** As requested, check MFEM implementation (e.g., `fem/integ.hpp`) to see base classes.
      * **Custom Class Definition (in the test file):**
        ```cpp
        // Test class that scales a mass matrix.
        class MyTestIntegrator : public mfem::BilinearFormIntegrator {
        private:
            mfem::ConstantCoefficient scale;
        public:
            explicit MyTestIntegrator(double s) : scale(s) {}
            void AssembleElementMatrix(const mfem::FiniteElement &el,
                                       mfem::ElementTransformation &Tr,
                                       mfem::DenseMatrix &elmat) override {
                mfem::MassIntegrator massInteg;
                massInteg.AssembleElementMatrix(el, Tr, elmat);
                elmat *= scale.constant;
            }
        };
        ```
      * **Test Logic:**
          * [ ] Setup a `FemProblem` with `PhysicsElectrostatics` (as in the test above).
          * [ ] Add the custom integrator: `auto customInteg = std::make_unique<MyTestIntegrator>(123.0);`
          * [ ] `physics->addBilinearIntegrator(std::move(customInteg));`
          * [ ] Add a source term, setup, and solve. Store the result `x_custom`.
          * [ ] **Verification:** Setup a *second* `FemProblem` (`problem2`).
          * [ ] Add a *standard* `mfem::MassIntegrator` with a coefficient of 123.0.
          * [ ] `auto *coeff = new mfem::ConstantCoefficient(123.0);`
          * [ ] `auto massInteg = std::make_unique<mfem::MassIntegrator>(*coeff);`
          * [ ] `physics2->addBilinearIntegrator(std::move(massInteg));`
          * [ ] Add the *same* source term, setup, and solve. Store the result `x_standard`.
          * [ ] **Assertion:** `ASSERT(mfem::Distance(x_custom, x_standard) < 1.e-12)`.
      * This test fails until Step 2.3 is complete.

#### Step 2.3: Implementation

  * **Checklist:**
      * [ ] Create and implement `src/hpcfem/integrator_info.hpp` as planned.
      * [ ] Modify `src/hpcfem/physics_interface.hpp` (add vectors and `add...` methods) as planned.
      * [ ] Implement `src/hpcfem/physics_interface.cpp` (`add...` methods) as planned.
      * [ ] Modify `src/hpcfem/physics_electrostatics.hpp / .cpp` (remove `set...`, coefficients, update `setupOn`) as planned.
      * [ ] Modify `src/hpcfem/physics_thermal.hpp / .cpp` (remove `set...`, coefficients, update `setupOn`) as planned.

#### Step 2.4: Testing, Refactoring, and Documentation

1.  **Testing:**
      * **Checklist:**
          * [ ] Run `tests/test_problem_abstraction.cpp`. It must pass.
          * [ ] Run `tests/test_custom_integrator.cpp`. It must pass.
2.  **Refactoring (Guideline 7):**
      * **Checklist:**
          * [ ] **Completely rewrite `example/ex0`, `ex1`, and `ex3`.**
          * [ ] **Hint:** These files are now the *user's* responsibility. They must manually create all `mfem::Coefficient`s (for epsilon, source, etc.), create all `mfem::...Integrator`s, and add them to the `physics` module using the new `add...Integrator` methods *before* calling `problem->addPhysics()`.
          * [ ] Verify all examples compile and run correctly.
3.  **Documentation (Guideline 5, 6, 22):**
      * **Checklist:**
          * [ ] Add new `struct` and method names to `docs/hpcfem-doc/naming_registry.md`.
          * [ ] Add Doxygen comments for all new `struct`s and `add...` methods.
          * [ ] Create `docs/mfem-basic-usage/user_defined_physics.md` explaining the new workflow (create fespace -\> create problem -\> create physics -\> create coeffs -\> create integrators -\> add integrators -\> add physics -\> setup -\> solve).
          * [ ] Update `docs/hpcfem-doc/custom_fem_interfaces.md` to mark Part 2 as complete.

-----

### Part 3: Parallel and Stress Testing (Guideline 10, 25)

**Goal:** Ensure the new, flexible design works correctly with MPI and scales to large problems.

#### Step 3.1: Write Failing Test (Parallel)

1.  **Create New Test: `tests/test_parallel_custom_integrator.cpp`**
      * **Checklist:**
          * [ ] Copy `tests/test_custom_integrator.cpp`.
          * [ ] **Hint:** This test *must* be runnable with `mpirun -np N ...`.
          * [ ] The test *must* initialize MPI (`mfem::Mpi::Init();`).
          * [ ] The test must use `mfem::ParMesh`, `mfem::ParFiniteElementSpace`, and `hpcfem::FemProblem`.
          * [ ] The test logic (comparing `MyTestIntegrator` to a standard `MassIntegrator`) remains identical.
          * [ ] **Assertion:** The parallel solution comparison `mfem::Distance(x_custom, x_standard)` must be `< 1.e-12` when run on 2, 4, or more processors.
      * This test will fail if any assumptions made in the refactoring were not parallel-safe.

#### Step 3.2: Implementation (Bug Fixes)

  * **Checklist:**
      * [ ] Run the `test_parallel_custom_integrator` with `mpirun -np 4`.
      * [ ] **Implementation:** No new code is expected. This step is purely for *fixing* any parallel bugs exposed by the test. The new design should be inherently parallel-safe as it relies on MFEM's parallel forms (`ParBilinearForm`, etc.), and the integrators are assembled element-wise.

#### Step 3.3: Stress Test

1.  **Create New Benchmark: `benchmark/poisson_scaling_femproblem/main.cpp`**
      * **Checklist:**
          * [ ] Create a new benchmark executable (add to `benchmark/CMakeLists.txt`).
          * [ ] This benchmark will mimic `benchmark/poisson_scaling`, but will use the newly refactored `FemProblem` and `PhysicsElectrostatics`.
          * [ ] The `main` function should contain a loop over `ref_levels`.
          * [ ] **Inside the loop:**
              * [ ] Create `pmesh` and refine it.
              * [ ] Create `fespace` (e.g., H1).
              * [ ] Create `problem` (FemProblem).
              * [ ] Create `physics` (PhysicsElectrostatics).
              * [ ] Create `coeff` and `integ` (e.g., `DiffusionIntegrator`).
              * [ ] `physics->addBilinearIntegrator(...)`.
              * [ ] (Add source term integrator).
              * [ ] `problem->addPhysics(physics.get());`
              * [ ] `problem->setup();`
              * [ ] `problem->solve();`
          * [ ] Run this benchmark with MPI (`mpirun`) and increase `ref_levels` until the DoFs exceed 1e7 (Guideline 10).
      * **Goal:** The test passes if it completes all refinement levels without crashing or producing errors.

#### Step 3.4: Final Documentation and Merge

  * **Checklist:**
      * [ ] (Guideline 7) Ensure *all* tests, examples, and benchmarks compile and run.
      * [ ] (Guideline 6) Ensure all new code has clean Doxygen comments.
      * [ ] (Guideline 22) Finalize `docs/hpcfem-doc/custom_fem_interfaces.md`.
      * [ ] (Guideline 8) Search the branch for any `TODO`s added and resolve them.
      * [ ] (Guideline 9) Merge the `feature/custom_fem_interfaces` branch into `dev` and delete the feature branch.