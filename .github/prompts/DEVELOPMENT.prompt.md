---
mode: agent
---

### **TDD Development Plan: General Waveguide Simulation Tool**

**Objective:** To create a general, robust, and parallel-capable tool for electromagnetic waveguide simulation within the `hpcfem` project. This tool will be developed in two main phases, both strictly following TDD.

1.  **Phase 1: 2D Eigenmode Solver:** Solves for the modal propagation constants (e.g., $k_c$, $\gamma$) and 2D vector field profiles (e.g., $TE_{10}$) on a waveguide cross-section. This is the foundation for the 3D solver's ports.
2.  **Phase 2: 3D Frequency-Domain (S-Parameter) Solver:** Solves the full 3D vector Helmholtz equation for a driven problem to calculate S-parameters.

---

### **Phase 1: 2D Vector Eigenmode Solver**

**Files to Create:**
* `src/hpcfem/physics/physics_waveguide_2d_eigen.hpp`
* `src/hpcfem/physics/physics_waveguide_2d_eigen.cpp`
* `tests/test_physics_waveguide_2d_eigen.cpp`

**Documentation:**
* Add `hpcfem::PhysicsWaveguide2DEigen` to `docs/hpcfem-doc/naming_registry.md`.
* Create `docs/hpcfem-doc/physics_waveguide_2d_eigen.md` explaining the physics and usage.

#### **TDD Step-by-Step Construction:**

**Step 1.1: Stub Class and Basic Test**
* **Test:** `test_class_instantiation`
    * Write a test in `test_physics_waveguide_2d_eigen.cpp` that includes the new header.
    * The test will instantiate `hpcfem::PhysicsWaveguide2DEigen`.
    * Assert that the pointer is not null.
* **Implementation:**
    * Create the stub `PhysicsWaveguide2DEigen` class in its `.hpp`/`.cpp` files.
    * It must inherit from the `hpcfem::PhysicsInterface`.
    * Implement empty virtual functions (e.g., `setup()`, `assemble()`, `solve()`).

**Step 1.2: Test 2D $H(curl)$ Space and PEC Boundary**
* **Test:** `test_nd_space_and_pec_boundary`
    * Load a 2D mesh (e.g., a simple 4-element rectangle).
    * Define a boundary attribute for the outer wall (e.g., attribute 1).
    * The test will:
        1.  Initialize the `PhysicsWaveguide2DEigen` and call its `setup()` method, passing in the mesh and the PEC boundary attribute.
        2.  Assert that the internal `FiniteElementSpace` is *not* null.
        3.  Assert that the space is of type `mfem::FiniteElementCollection::NEDELEC_FE`.
        4.  Assert that the `FiniteElementSpace::GetEssentialTrueDofs()` list is non-empty and contains the correct number of edge DoFs corresponding to the outer boundary.
* **Implementation:**
    * In `PhysicsWaveguide2DEigen::setup()`, create and store a `mfem::ParFiniteElementSpace`.
    * Use an `mfem::ND_FECollection` (Nedelec $H(curl)$ elements).
    * Store the PEC boundary attribute list (`mfem::Array<int>`).
    * Call `fespace->GetEssentialTrueDofs(pecBdr, essTDoFs)`.

**Step 1.3: Test Bilinear Form Assembly (Lossless Case)**
* **Test:** `test_eigen_matrix_assembly_lossless`
    * Use a known 2x1 mesh.
    * Call the `assemble()` method.
    * The test must get the assembled operators (A and B).
    * Assert that both `mfem::HypreParMatrix` pointers (for A and B) are not null.
    * Assert that their dimensions are correct (equal to `fespace->GetTrueVSize()`).
    * **Physics:** We are solving $\nabla_t \times (\mu_r^{-1} \nabla_t \times \vec{E}_t) = k_c^2 \epsilon_r \vec{E}_t$.
* **Implementation:**
    * `assemble()` will create two `mfem::ParBilinearForm`s: `a` and `b`.
    * `a->AddDomainIntegrator(new mfem::CurlCurlIntegrator(muInvCoeff));` (This is $[S_t]$).
    * `b->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(epsilonCoeff));` (This is $[T_t]$).
    * `a->Assemble()` and `b->Assemble()`.
    * Finalize the operators using `a->ParallelAssemble()` and `b->ParallelAssemble()`, and apply essential BCs.

**Step 1.4: Test Eigenvalue Solution (Lossless $TE_{10}$)**
* **Test:** `test_te10_cutoff_frequency_serial_and_parallel`
    * Load a mesh of a standard rectangular waveguide (e.g., WR-90, width `a = 22.86e-3` m).
    * Run the full `setup()`, `assemble()`, and `solve()`.
    * **Check:** The analytical cutoff wavenumber is $k_c = \pi / a$. The eigenvalue $\lambda$ from the solver is $k_c^2$.
    * The test must:
        1.  Call `solve()` and get the lowest eigenvalue $\lambda_0$.
        2.  Calculate the expected $\lambda_{exp} = (\pi / 0.02286)^2$.
        3.  Assert that `abs(lambda_0 - lambda_exp) < 1e-3`.
    * This test *must* be run in serial (`mpirun -n 1 ...`) and parallel (`mpirun -n 4 ...`) and pass both.
* **Implementation:**
    * You will need an eigenvalue solver. Add `src/hpcfem/solvers/solver_eigen_slepc.hpp` (or similar) that wraps the SLEPc library.
    * `PhysicsWaveguide2DEigen::solve()` will call this solver to solve $A\{e\} = \lambda B\{e\}$.
    * The solver should return the computed eigenvalues.

**Step 1.5: Test Lossy Media (Complex Eigenvalue)**
* **Test:** `test_lossy_waveguide_complex_kc`
    * Use the same WR-90 mesh.
    * In the test, define a *complex* permittivity $\epsilon_r = \epsilon' - j\epsilon''$ (where $\epsilon'' > 0$).
    * Pass this as a `mfem::ComplexConstantCoefficient` to the physics class.
    * Run the full `setup()`, `assemble()`, and `solve()`.
    * **Check:** The eigenvalue $\lambda = k_c^2$ must now be complex.
    * Assert that the imaginary part of the computed eigenvalue $\lambda_0$ is non-zero.
* **Implementation:**
    * `PhysicsWaveguide2DEigen` must be updated to handle `mfem::ComplexOperator`.
    * The `assemble()` method must create complex-valued `ParBilinearForm`s.
    * The `SolverEigenSlepc` must be configured for complex arithmetic.

---

### **Phase 2: 3D Frequency-Domain (S-Parameter) Solver**

**Files to Create:**
* `src/hpcfem/physics/waveguide_port.hpp` (a helper class)
* `src/hpcfem/physics/waveguide_port.cpp`
* `src/hpcfem/physics/physics_waveguide_3d_freq.hpp`
* `src/hpcfem/physics/physics_waveguide_3d_freq.cpp`
* `src/hpcfem/physics/port_excitation_integrator.hpp` (custom integrator)
* `tests/test_physics_waveguide_3d_freq.cpp`

**Documentation:**
* Add new classes to `naming_registry.md`.
* Create `docs/hpcfem-doc/physics_waveguide_3d_freq.md` explaining S-parameter simulation.

#### **TDD Step-by-Step Construction:**

**Step 2.1: Test Port Helper Class**
* **Test:** `test_waveguide_port_initialization` (in `test_physics_waveguide_3d_freq.cpp` or a new file)
    * This test must first run the **Phase 1 Solver** on a 2D port mesh to get the $TE_{10}$ mode (field and $\beta$).
    * It then instantiates `hpcfem::WaveguidePort` with this modal solution.
    * Assert that the port object correctly stores the modal field and propagation constant.
* **Implementation:**
    * Create `WaveguidePort` class. It stores the 2D `FiniteElementSpace`, the modal `GridFunction` $\vec{e}_1(x,y)$, and the complex propagation constant $\gamma_1 = \alpha + j\beta_1$.

**Step 2.2: Test 3D Matrix Assembly ($[S] - k_0^2 [M]$)**
* **Test:** `test_3d_system_matrix_assembly`
    * Load a 3D mesh (e.g., `testmesh_cube.mesh`).
    * Set up `PhysicsWaveguide3DFreq` with a frequency $k_0$.
    * Call `assemble()`.
    * **Check:** Assert that the final system matrix $[A]$ is a non-null `mfem::HypreParMatrix` and is complex.
* **Implementation:**
    * `PhysicsWaveguide3DFreq::assemble()` will create the complex system matrix:
        * `s = new mfem::ParBilinearForm(fespace, mfem::ComplexOperator::HERMITIAN);`
        * `s->AddDomainIntegrator(new mfem::CurlCurlIntegrator(muInvCoeff));`
        * `m = new mfem::ParBilinearForm(fespace, mfem::ComplexOperator::HERMITIAN);`
        * `m->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(epsilonCoeff));`
        * `s->Assemble(); m->Assemble();`
        * `A = s + (-k0*k0)*m;` (Symbolic assembly)
        * `A->Finalize();`

**Step 2.3: Test Port Excitation Vector $\{b\}$**
* **Test:** `test_port_excitation_rhs`
    * Load a 3D cube mesh. Designate one face (e.g., attribute 1) as Port 1.
    * Create a `WaveguidePort` object (from 2.1) and assign it to attribute 1.
    * Call `assembleRHS()` (or equivalent).
    * **Check:**
        1.  Get the assembled `mfem::HypreParVector` $\{b\}$.
        2.  Assert that $\{b\}$ is non-zero.
        3.  Manually check that the non-zero entries correspond *only* to DoFs on the port surface.
* **Implementation:**
    * `PhysicsWaveguide3DFreq::assembleRHS()` will create a `mfem::ParLinearForm`.
    * You must create a custom integrator: `PortExcitationIntegrator : mfem::BoundaryLFIntegrator`.
    * This integrator's `AssembleRHSElementVect` will take the $TE_{10}$ modal field $\vec{e}_1$ (as a `VectorCoefficient`) and compute the integral: $\int_{S_{p1}} \vec{W}_i \cdot (2 j \beta_1 Y_1 \vec{e}_1) dS$.
    * Add this integrator to the `ParLinearForm`: `b->AddBoundaryIntegrator(new PortExcitationIntegrator(...), port1Attr);`
    * `b->Assemble();`

**Step 2.4: Test Port Absorbing Boundary (on $[A]$)**
* **Test:** `test_port_absorption_lhs`
    * Use the same setup as 2.3.
    * Call `assemble()` (which now must also assemble port BCs).
    * **Check:** The assembled matrix $[A]$ must be different from the one in test 2.2.
* **Implementation:**
    * The port BC (first-order) adds an impedance term to the $[A]$ matrix:
    * $\int_{S_p} \vec{W}_t \cdot (Y_1 \vec{E}_t) dS$ (where $Y_1$ is mode 1 impedance).
    * `PhysicsWaveguide3DFreq::assemble()` must now *also* add this boundary term to the `s` or `m` form:
    * `a->AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(portImpedanceCoeff), portAttr);`
    * *Note:* A true modal absorbing boundary is more complex, but this $TE_{10}$-matched impedance is a valid TDD starting point.

**Step 2.5: Test S-Parameter Calculation (S11/S21)**
* **Test:** `test_s_params_straight_waveguide_serial_and_parallel`
    * **Problem:** A 3D mesh of a straight waveguide segment.
    * **Setup:**
        1.  Port 1 (e.g., attr 1, $z=0$) - Excite + Absorb.
        2.  Port 2 (e.g., attr 2, $z=L$) - Absorb Only.
        3.  PEC walls (e.g., attr 3).
    * Run the full `setup()`, `assemble()`, `solve()`.
    * Call a new `calculateSParameters()` post-processing function.
    * **Check (Lossless Case):**
        1.  Assert `abs(S11) < 1e-3` (low reflection).
        2.  Assert `abs(abs(S21) - 1.0) < 1e-3` (low loss).
        3.  Assert `abs(S11_mag_db) < -60.0`.
    * This test *must* pass in serial and parallel.
* **Implementation:**
    * `PhysicsWaveguide3DFreq::solve()` will use a linear solver (e.g., `SolverHypreAmg` or PETSc) to solve $A\{e\} = \{b\}$ for the complex field $\{e\}$.
    * Create `PhysicsWaveguide3DFreq::calculateSParameters()`.
    * This function needs the total field solution `mfem::ParComplexGridFunction e_total`.
    * **S11:**
        1.  Project $e_{\text{total}}$ onto the $TE_{10}$ mode $\vec{e}_1$ at Port 1:
            $C_1 = \text{project}(\vec{E}_{\text{total}}, \vec{e}_1, \text{Port 1})$.
        2.  $S_{11} = C_1 - 1.0$ (since incident amplitude $a_1=1.0$).
    * **S21:**
        1.  Project $e_{\text{total}}$ onto the $TE_{10}$ mode $\vec{e}_1$ at Port 2:
            $C_2 = \text{project}(\vec{E}_{\text{total}}, \vec{e}_1, \text{Port 2})$.
        2.  $S_{21} = C_2$.
    * The projection integral is: $C = (\int_{S_p} \vec{E}_{\text{total}} \cdot \vec{e}_1^* dS) / (\int_{S_p} \vec{e}_1 \cdot \vec{e}_1^* dS)$. You will need `mfem::BilinearForm` and `mfem::LinearForm` with `VectorFEMassIntegrator` on the *boundary* to compute these.

**Step 2.6: Final Refactor and Stress Test**
* **Test:** `test_discontinuity_power_conservation`
    * **Problem:** A 3D waveguide with a dielectric block or iris (a real device).
    * **Check (Lossless):** Assert that $|S_{11}|^2 + |S_{21}|^2$ is within tolerance of $1.0$.
* **Stress Test:**
    * Rerun `test_s_params_discontinuity` on a mesh with >1M DoFs.
    * Rerun in parallel on >16 cores.
    * The test must still pass, and memory usage should be reasonable.
* **Refactor:** Clean up all code, ensure all Doxygen comments are complete, and clear all `TODO`s before merging to `dev`.