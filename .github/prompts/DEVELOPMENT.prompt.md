---
mode: agent
---

# HPC-FEM Development Guidelines and Roadmap

Here is a comprehensive TDD roadmap to guide an AI coder in polishing the `hpc-fem-playground` library. This plan focuses on (1) refactoring the `src` directory structure and (2) implementing a generic interface for user-defined finite element spaces and integrators, all while adhering strictly to your provided developing guidelines.

-----

## TDD Roadmap: `hpc-fem-playground` Refinement

This roadmap is divided into two main phases. Each step must follow the TDD workflow (Guideline 10: PLAN -\> TEST -\> IMPLEMENT -\> TEST -\> REFACTOR -\> DOCUMENT). All code must adhere to the **DEVELOPING GUIDELINES**.

### Part 1: Clarify `src` File Structure

**Goal:** Refactor the `src/hpcfem` directory into a more organized, scalable structure by categorizing files. The "test" for this refactor is that the library compiles and all existing tests in the `tests/` directory pass after each step.

**New Target Structure:**

```
src/
├── CMakeLists.txt
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
        ├── solver_block_gauss_seidel.hpp
        └── solver_block_gauss_seidel.cpp
```

-----

#### Step 1.1: Refactor Solvers

  * **PLAN:** Move all solver-related files (`solver_*.hpp`, `solver_*.cpp`) into a new `src/hpcfem/solvers/` directory.
  * **IMPLEMENT (Refactor):**
    1.  Create directory `src/hpcfem/solvers/`.
    2.  Move `src/hpcfem/solver_hypre_amg.hpp`, `src/hpcfem/solver_hypre_amg.cpp`, `src/hpcfem/solver_block_gauss_seidel.hpp`, and `src/hpcfem/solver_block_gauss_seidel.cpp` into the new directory.
    3.  Move `src/hpcfem/solver_interface.hpp` into `src/hpcfem/core/` (create this directory first).
  * **IMPLEMENT (Fix):**
    1.  Update `src/CMakeLists.txt` to reflect the new file paths for the `hpcfem` library target.
        ```cmake
        # Hint: Change file paths from 'hpcfem/solver_...' to 'hpcfem/solvers/solver_...'
        add_library(hpcfem
            hpcfem/fem_problem.cpp
            ...
            hpcfem/solvers/solver_hypre_amg.cpp
            hpcfem/solvers/solver_block_gauss_seidel.cpp
        )
        ```
    2.  Update all `#include` directives in all files (e.g., `fem_problem.cpp`, `physics_joule_heating.cpp`) to use the new paths. Per Guideline 18, paths should be relative to `src`.
        ```cpp
        // Old: #include "hpcfem/solver_interface.hpp"
        // New: #include "hpcfem/core/solver_interface.hpp"

        // Old: #include "hpcfem/solver_hypre_amg.hpp"
        // New: #include "hpcfem/solvers/solver_hypre_amg.hpp"
        ```
  * **TEST:**
    1.  Re-run CMake configuration.
    2.  Compile the entire project (library, tests, examples).
    3.  Run all executables in the `tests/` directory.
  * **CHECKLIST:**
      * [ ] `src/hpcfem/solvers/` directory created.
      * [ ] `src/hpcfem/core/` directory created.
      * [ ] Solver files moved to `solvers/`.
      * [ ] `solver_interface.hpp` moved to `core/`.
      * [ ] `src/CMakeLists.txt` updated with new paths.
      * [ ] All `#include` paths updated project-wide.
      * [ ] Project compiles successfully.
      * [ ] All existing tests pass.

-----

#### Step 1.2: Refactor Physics

  * **PLAN:** Move all physics-related files (`physics_*.hpp`, `physics_*.cpp`) into a new `src/hpcfem/physics/` directory.
  * **IMPLEMENT (Refactor):**
    1.  Create directory `src/hpcfem/physics/`.
    2.  Move `src/hpcfem/physics_electrostatics.hpp`, `src/hpcfem/physics_electrostatics.cpp`, `src/hpcfem/physics_thermal.hpp`, `src/hpcfem/physics_thermal.cpp`, `src/hpcfem/physics_joule_heating.hpp`, and `src/hpcfem/physics_joule_heating.cpp` into the new directory.
    3.  Move `src/hpcfem/physics_interface.hpp` into `src/hpcfem/core/`.
  * **IMPLEMENT (Fix):**
    1.  Update `src/CMakeLists.txt` to reflect the new file paths.
    2.  Update all `#include` directives (e.g., in `fem_problem.cpp`, `fem_problem.hpp`, and example files).
        ```cpp
        // Old: #include "hpcfem/physics_interface.hpp"
        // New: #include "hpcfem/core/physics_interface.hpp"

        // Old: #include "hpcfem/physics_electrostatics.hpp"
        // New: #include "hpcfem/physics/physics_electrostatics.hpp"
        ```
  * **TEST:**
    1.  Compile the entire project.
    2.  Run all executables in the `tests/` directory.
  * **CHECKLIST:**
      * [ ] `src/hpcfem/physics/` directory created.
      * [ ] Physics files moved to `physics/`.
      * [ ] `physics_interface.hpp` moved to `core/`.
      * [ ] `src/CMakeLists.txt` updated.
      * [ ] All `#include` paths updated project-wide.
      * [ ] Project compiles successfully.
      * [ ] All existing tests pass.

-----

#### Step 1.3: Relocate Core Files

  * **PLAN:** Move the remaining `fem_problem` files into the `src/hpcfem/core/` directory.
  * **IMPLEMENT (Refactor):**
    1.  Move `src/hpcfem/fem_problem.hpp` and `src/hpcfem/fem_problem.cpp` into `src/hpcfem/core/`.
  * **IMPLEMENT (Fix):**
    1.  Update `src/CMakeLists.txt`.
    2.  Update all `#include` directives (mainly in `tests/` and `example/` files).
        ```cpp
        // Old: #include "hpcfem/fem_problem.hpp"
        // New: #include "hpcfem/core/fem_problem.hpp"
        ```
  * **TEST:**
    1.  Compile the entire project.
    2.  Run all executables in the `tests/` and `example/` directories.
  * **DOCUMENT:**
    1.  Create a new document `docs/hpcfem-doc/file_structure.md`.
    2.  Document the new file structure.
  * **CHECKLIST:**
      * [ ] `fem_problem` files moved to `core/`.
      * [ ] `src/CMakeLists.txt` updated.
      * [ ] All `#include` paths updated project-wide.
      * [ ] Project compiles successfully.
      * [ ] All existing tests and examples run successfully.
      * [ ] `docs/hpcfem-doc/file_structure.md` created and populated.

-----

-----

### Part 2: Add Interface for User-Defined FE Space and Integrators

**Goal:** Create a generic physics abstraction that allows users to define a problem by supplying a custom `FiniteElementSpace` and lists of `Integrator` wrappers. This allows building new physics problems without writing a new `Physics` class each time.

-----

#### Step 2.1: `FiniteElementSpace` Abstraction

  * **PLAN:** Create a class `hpcfem::FiniteElementSpace` to encapsulate the creation and management of `mfem::FiniteElementSpace` and its associated `mfem::FiniteElementCollection`.
  * **TEST:** Create `tests/test_fe_space.cpp`.
      * `test_h1_space_creation_serial`: Create a serial `mfem::Mesh`, then create `hpcfem::FiniteElementSpace(mesh, 1, 3, FiniteElementSpace::Type::H1)`. Verify `getMfemSpace()->GetOrder()` is 1 and `getMfemSpace()->GetVDim()` is 1.
      * `test_nd_space_creation_serial`: Create a serial `mfem::Mesh`, then create `hpcfem::FiniteElementSpace(mesh, 1, 3, FiniteElementSpace::Type::ND_VECTOR)`. Verify `getMfemSpace()->GetOrder()` is 1.
      * `test_parallel_space_creation`: (Compile-time guarded with `#ifdef MFEM_USE_MPI`) Create a `mfem::ParMesh`, then create an `hpcfem::FiniteElementSpace`. Verify `getMfemSpace()` returns a valid `mfem::ParFiniteElementSpace`.
  * **IMPLEMENT:** Create `src/hpcfem/core/finite_element_space.hpp` and `.cpp`.
      * `finite_element_space.hpp`:
        ```cpp
        #ifndef HPCFEM_FINITE_ELEMENT_SPACE_HPP
        #define HPCFEM_FINITE_ELEMENT_SPACE_HPP

        #include "mfem.hpp"

        namespace hpcfem
        {

        class FiniteElementSpace
        {
        public:
            enum class Type { H1, L2, ND_VECTOR, RT_VECTOR };

        #ifdef MFEM_USE_MPI
            FiniteElementSpace(mfem::ParMesh* pmesh, int order, int dim = 1, Type type = Type::H1);
            mfem::ParFiniteElementSpace* getMfemSpace() const;
        #else
            FiniteElementSpace(mfem::Mesh* mesh, int order, int dim = 1, Type type = Type::H1);
            mfem::FiniteElementSpace* getMfemSpace() const;
        #endif

            ~FiniteElementSpace();

            mfem::FiniteElementCollection* getFeCollection() const;

        private:
            mfem::FiniteElementCollection* fec_ = nullptr; // Owned
        #ifdef MFEM_USE_MPI
            mfem::ParFiniteElementSpace* fes_ = nullptr; // Owned
            mfem::ParMesh* pmesh_ = nullptr; // Not owned
        #else
            mfem::FiniteElementSpace* fes_ = nullptr; // Owned
            mfem::Mesh* mesh_ = nullptr; // Not owned
        #endif
            bool isOwner_ = true; // Manages fec_ and fes_
        };

        } // namespace hpcfem
        #endif // HPCFEM_FINITE_ELEMENT_SPACE_HPP
        ```
      * `finite_element_space.cpp`: Implement the constructor to `new` the appropriate `mfem::FiniteElementCollection` (e.g., `mfem::H1_FECollection`) and `mfem::FiniteElementSpace` based on `type`, `order`, and `dim`. Implement the destructor to `delete fec_` and `fes_`.
  * **TEST (Fix):** Compile and run `test_fe_space`. Debug until all tests pass.
  * **DOCUMENT:** Add Doxygen comments. Register `FiniteElementSpace` and `FiniteElementSpace::Type` in `docs/hpcfem-doc/naming_registry.md`.
  * **CHECKLIST:**
      * [ ] `tests/test_fe_space.cpp` created with all test cases.
      * [ ] `src/hpcfem/core/finite_element_space.hpp` created.
      * [ ] `src/hpcfem/core/finite_element_space.cpp` created.
      * [ ] `src/CMakeLists.txt` and `tests/CMakeLists.txt` updated for new files.
      * [ ] All tests in `test_fe_space` pass.
      * [ ] Documentation and naming registry updated.

-----

#### Step 2.2: `Integrator` Interfaces

  * **PLAN:** Define abstract base classes for different types of integrators (bilinear form, linear form, etc.). These will be simple wrappers.
  * **TEST:** Create `tests/test_integrator_interfaces.cpp`.
      * `test_mock_integrators`: Create mock classes inheriting from each interface (e.g., `class MockBilinear : public BilinearIntegratorInterface {}`) and instantiate them to ensure the interfaces are correctly defined.
  * **IMPLEMENT:** Create `src/hpcfem/core/integrator_interface.hpp`.
    ```cpp
    #ifndef HPCFEM_INTEGRATOR_INTERFACE_HPP
    #define HPCFEM_INTEGRATOR_INTERFACE_HPP

    #include "mfem.hpp"

    namespace hpcfem
    {

    /** @brief Base for all bilinear form integrators (domain and boundary) */
    class BilinearIntegratorInterface
    {
    public:
        virtual ~BilinearIntegratorInterface() = default;
        /** @brief Adds the specific MFEM integrator to the form.
         * CODER: Check mfem/fem/bilininteg.hpp */
        virtual void addToBilinearForm(mfem::BilinearForm& form) const = 0;
    };

    /** @brief Base for all linear form domain integrators */
    class LinearDomainIntegratorInterface
    {
    public:
        virtual ~LinearDomainIntegratorInterface() = default;
        /** @brief Adds the specific MFEM integrator to the form.
         * CODER: Check mfem/fem/lininteg.hpp */
        virtual void addToLinearForm(mfem::LinearForm& form) const = 0;
    };

    /** @brief Base for all linear form boundary integrators */
    class LinearBoundaryIntegratorInterface
    {
    public:
        virtual ~LinearBoundaryIntegratorInterface() = default;
        /** @brief Adds the specific MFEM integrator to the form.
         * CODER: Check mfem/fem/lininteg.hpp */
        virtual void addToLinearForm(mfem::LinearForm& form, const mfem::Array<int>& markers) const = 0;
    };

    } // namespace hpcfem
    #endif // HPCFEM_INTEGRATOR_INTERFACE_HPP
    ```
  * **TEST (Fix):** Compile and run `test_integrator_interfaces`.
  * **DOCUMENT:** Add Doxygen comments. Register all interfaces in `naming_registry.md`.
  * **CHECKLIST:**
      * [ ] `tests/test_integrator_interfaces.cpp` created and passes.
      * [ ] `src/hpcfem/core/integrator_interface.hpp` created.
      * [ ] `src/CMakeLists.txt` and `tests/CMakeLists.txt` updated.
      * [ ] Documentation and naming registry updated.

-----

#### Step 2.3: Concrete `Integrator` Wrappers

  * **PLAN:** Implement concrete wrappers for a few common MFEM integrators (e.g., `DiffusionIntegrator`, `DomainLFIntegrator`). **Coder must check `share/mfem/` examples** (like `ex1.cpp`) to see how these are used.
  * **TEST:** Create `tests/test_integrator_wrappers.cpp`.
      * `test_diffusion_integrator`: Create `mfem::ConstantCoefficient`, pass it to `hpcfem::DiffusionIntegrator`. Create a `mfem::BilinearForm` and call `addToBilinearForm()`. Check that the form is modified.
      * `test_domain_lf_integrator`: Create `mfem::ConstantCoefficient`, pass it to `hpcfem::DomainLFIntegrator`. Create a `mfem::LinearForm` and call `addToLinearForm()`. Check that the form is modified.
  * **IMPLEMENT:** Create `src/hpcfem/integrators/` directory.
      * `src/hpcfem/integrators/diffusion_integrator.hpp`:
        ```cpp
        // ... includes ...
        #include "hpcfem/core/integrator_interface.hpp"

        namespace hpcfem {
        class DiffusionIntegrator : public BilinearIntegratorInterface {
        public:
            /** @param coeff Not owned, must persist for life of this object */
            explicit DiffusionIntegrator(mfem::Coefficient& coeff) : coeff_(&coeff) {}
            void addToBilinearForm(mfem::BilinearForm& form) const override {
                form.AddDomainIntegrator(new mfem::DiffusionIntegrator(*coeff_));
            }
        private:
            mfem::Coefficient* coeff_;
        };
        } // namespace hpcfem
        ```
      * `src/hpcfem/integrators/domain_lf_integrator.hpp`:
        ```cpp
        // ... includes ...
        #include "hpcfem/core/integrator_interface.hpp"

        namespace hpcfem {
        class DomainLFIntegrator : public LinearDomainIntegratorInterface {
        public:
            /** @param coeff Not owned, must persist for life of this object */
            explicit DomainLFIntegrator(mfem::Coefficient& coeff) : coeff_(&coeff) {}
            void addToLinearForm(mfem::LinearForm& form) const override {
                form.AddDomainIntegrator(new mfem::DomainLFIntegrator(*coeff_));
            }
        private:
            mfem::Coefficient* coeff_;
        };
        } // namespace hpcfem
        ```
  * **TEST (Fix):** Compile and run `test_integrator_wrappers`.
  * **DOCUMENT:** Doxygen, `naming_registry.md`.
  * **CHECKLIST:**
      * [ ] `src/hpcfem/integrators/` directory created.
      * [ ] `tests/test_integrator_wrappers.cpp` created and passes.
      * [ ] `diffusion_integrator.hpp` created and implemented.
      * [ ] `domain_lf_integrator.hpp` created and implemented.
      * [ ] `src/CMakeLists.txt` updated (no new `.cpp` files, but good to add dir).
      * [ ] `tests/CMakeLists.txt` updated.
      * [ ] Documentation and naming registry updated.

-----

#### Step 2.4: `GenericPhysics` Implementation

  * **PLAN:** Create `hpcfem::GenericPhysics`, a concrete class inheriting from `PhysicsInterface`, that assembles a system based on lists of integrator wrappers.
  * **TEST:** Create `tests/test_generic_physics.cpp`.
      * `test_poisson_assembly_serial`:
        1.  Create serial mesh, `hpcfem::FiniteElementSpace` (H1, order 1).
        2.  Create `mfem::ConstantCoefficient k(1.0)` and `mfem::ConstantCoefficient f(1.0)`.
        3.  Create `hpcfem::DiffusionIntegrator diff(k)` and `hpcfem::DomainLFIntegrator source(f)`.
        4.  Create `std::vector`s holding pointers to these integrators.
        5.  Create `hpcfem::GenericPhysics` with the FE space and integrator lists.
        6.  Call `physics.assemble(A, b)`.
        7.  Check `A->Height()` and `b.Size()` match `fes->getMfemSpace()->GetTrueVSize()`.
        8.  Check that `A->NumNonZeroElems()` \> 0.
  * **IMPLEMENT:** Create `src/hpcfem/physics/generic_physics.hpp` and `.cpp`.
      * `generic_physics.hpp`:
        ```cpp
        // ... includes ...
        #include "hpcfem/core/physics_interface.hpp"
        #include "hpcfem/core/finite_element_space.hpp"
        #include "hpcfem/core/integrator_interface.hpp"

        namespace hpcfem {
        class GenericPhysics : public PhysicsInterface {
        public:
            GenericPhysics(
                FiniteElementSpace* fes,
                const std::vector<BilinearIntegratorInterface*>& bilinearInteg,
                const std::vector<LinearDomainIntegratorInterface*>& linearDomainInteg,
                const std::vector<LinearBoundaryIntegratorInterface*>& linearBoundaryInteg
            );
            
            // Override assemble methods from PhysicsInterface
            // ... (serial and parallel)
            void assemble(mfem::SparseMatrix*& A, mfem::Vector*& b) override;
            // ... add MPI version ...
            
        private:
            FiniteElementSpace* fes_; // Not owned
            std::vector<BilinearIntegratorInterface*> bilinearIntegrators_; // Not owned
            std::vector<LinearDomainIntegratorInterface*> linearDomainIntegrators_; // Not owned
            std::vector<LinearBoundaryIntegratorInterface*> linearBoundaryIntegrators_; // Not owned
            
            // Boundary attributes for linear boundary integrators
            // TODO: This needs to be expanded for practical use
        };
        } // namespace hpcfem
        ```
      * `generic_physics.cpp`:
        ```cpp
        // ... In assemble() method ...
        mfem::FiniteElementSpace* mfemFes = fes_->getMfemSpace();

        // 1. Create forms
        mfem::BilinearForm a(mfemFes);
        mfem::LinearForm b(mfemFes);

        // 2. Add Bilinear Integrators
        for (auto* integ : bilinearIntegrators_) {
            integ->addToBilinearForm(a);
        }

        // 3. Add Linear Domain Integrators
        for (auto* integ : linearDomainIntegrators_) {
            integ->addToLinearForm(b);
        }

        // 4. Add Linear Boundary Integrators
        // TODO: Handle boundary markers properly. This is a simplification.
        // mfem::Array<int> all_boundaries(pmesh_->bdr_attributes.Max());
        // all_boundaries = 1; 
        // for (auto* integ : linearBoundaryIntegrators_) {
        //     integ->addToLinearForm(b, all_boundaries);
        // }

        // 5. Assemble forms
        a.Assemble();
        b.Assemble();

        // 6. Apply Essential BCs
        // TODO: This interface needs to be extended to accept BCs
        mfem::Array<int> essTdofList; // Empty for now
        a.Finalize(essTdofList);

        // 7. Get system
        a.FormSystemMatrix(essTdofList, A); // A is output
        b.ParallelAssemble(b_out); // b_out is output
        ```
  * **TEST (Fix):** Compile and run `test_generic_physics`. Debug until pass.
  * **DOCUMENT:** Doxygen, `naming_registry.md`. Add `TODO` in `generic_physics.cpp` for boundary conditions.
  * **CHECKLIST:**
      * [ ] `tests/test_generic_physics.cpp` created and passes.
      * [ ] `src/hpcfem/physics/generic_physics.hpp` created.
      * [ ] `src/hpcfem/physics/generic_physics.cpp` created and implemented.
      * [ ] `src/CMakeLists.txt` and `tests/CMakeLists.txt` updated.
      * [ ] Assembly logic iterates over integrator lists.
      * [ ] Documentation and naming registry updated.

-----

#### Step 2.5: `FemProblem` Integration Test and Example

  * **PLAN:** Use the new `GenericPhysics` with the existing `FemProblem` class to solve a simple Poisson problem.
  * **TEST:** Create `tests/test_problem_generic_physics.cpp`.
      * `test_poisson_problem_end_to_end`:
        1.  Set up mesh, `hpcfem::FiniteElementSpace`, coefficients, and `hpcfem` integrators as in Step 2.4.
        2.  Create `hpcfem::GenericPhysics`.
        3.  Create a solver: `hpcfem::SolverHypreAmg solver;`
        4.  Create `hpcfem::FemProblem problem(mesh, &physics, &solver);`.
        5.  Call `problem.assemble()`.
        6.  Call `problem.solve()`.
        7.  Get solution `const mfem::Vector& x = problem.getSolution();`.
        8.  Check `x.Norml2()` is greater than 0 and not NaN.
  * **IMPLEMENT:**
    1.  No new library code. This step just uses the new components.
    2.  Create `example/ex7_generic_poisson.cpp` based on the test. Ensure it compiles and runs, saving a solution file.
  * **TEST (Fix):** Compile and run `test_problem_generic_physics` and `ex7_generic_poisson`.
  * **DOCUMENT:** Clean up Doxygen comments.
  * **CHECKLIST:**
      * [ ] `tests/test_problem_generic_physics.cpp` created and passes.
      * [ ] `example/ex7_generic_poisson.cpp` created.
      * [ ] `example/CMakeLists.txt` updated for `ex7`.
      * [ ] `ex7` compiles, runs, and produces a valid result.
      * [ ] All TODOs are documented.