---
mode: agent
---

# TDD Plan: General Waveguide Simulation Tool

**Project:** `hpc-fem-playground`
**Branch Name:** `feature/waveguide_toolkit`
**Goal:** To build a general, 3D, frequency-domain tool for waveguide simulation (S-Parameters) using MFEM's vector finite elements (Nedelec/edge elements).

This plan follows the TDD workflow (Rule \#10: PLANNING -\> TESTS -\> IMPLEMENTING -\> ...).

-----

## Step 1: The Core Physics - Waveguide Eigenvalue Problem

**1. Planning:**
Before we can solve a *driven* problem for S-parameters, we must be able to solve the *eigenvalue* problem. This will:

1.  Verify our use of vector finite elements (`ND_FECollection`) is correct.
2.  Prove we have eliminated the spurious mode (non-positive-definite matrix) problem, as discussed.
3.  Give us the "basis" (the modes) we will need later for port definitions.

We will solve the 2D eigenvalue problem for the transverse magnetic field $\mathbf{H}_t$ on a rectangular waveguide cross-section. The equation is $\nabla_t \times (\frac{1}{\epsilon_r} \nabla_t \times \mathbf{H}_t) = k_c^2 \mu_r \mathbf{H}_t$, where $k_c$ is the cutoff wavenumber (the eigenvalue).

**2. Test (WRITE THIS FIRST):**
Create a new file: `tests/test_physics_waveguide_eigen.cpp`

```cpp
// in tests/test_physics_waveguide_eigen.cpp
#include "gtest/gtest.h"
#include "mfem.hpp"
#include "hpcfem/physics/physics_waveguide_eigen.hpp"
#include "hpcfem/core/fem_problem.hpp"

#include <string>
#include <vector>
#include <cmath>

// Test case for Rule #10 and #20
TEST(PhysicsWaveguideEigen, test_rectangular_waveguide_cutoff) {
    // 1. Setup the problem
    constexpr int dim = 2;
    constexpr int order = 1;
    std::string meshFile = "../testdata/testmesh_rectangle_2d.mesh"; // TODO: Create this 2D mesh
    
    // Waveguide dimensions (e.g., WR-90)
    constexpr double a = 2.286e-2; // meters
    constexpr double b = 1.016e-2; // meters
    constexpr int numModes = 5;

    // 2. Create the FEM problem
    hpcfem::FemProblem problem(dim, order, meshFile);
    
    // 3. Create and set the physics
    hpcfem::PhysicsWaveguideEigen physics(problem.getMesh());
    problem.setPhysics(&physics);
    
    // 4. Solve the eigenvalue problem
    // TODO: The physics object will need a 'solve' method that returns eigenvalues
    std::vector<double> kcSq = physics.solve(numModes); // Solves for eigenvalues (k_c^2)

    // 5. Assertions: Check against analytical solution
    // k_c^2 = (m*pi/a)^2 + (n*pi/b)^2
    // We expect the first few modes to be TE10, TE20, TE01, TE30, TE11 (or similar)
    // Note: FEM solvers find them in ascending order, not by (m,n) index.
    
    constexpr double PI = 3.141592653589793;
    std::vector<double> kcSqAnalytical;
    kcSqAnalytical.push_back(std::pow(1*PI/a, 2) + std::pow(0*PI/b, 2)); // TE10
    kcSqAnalytical.push_back(std::pow(2*PI/a, 2) + std::pow(0*PI/b, 2)); // TE20
    kcSqAnalytical.push_back(std::pow(0*PI/a, 2) + std::pow(1*PI/b, 2)); // TE01
    kcSqAnalytical.push_back(std::pow(3*PI/a, 2) + std::pow(0*PI/b, 2)); // TE30
    kcSqAnalytical.push_back(std::pow(1*PI/a, 2) + std::pow(1*PI/b, 2)); // TE11
    
    std::sort(kcSqAnalytical.begin(), kcSqAnalytical.end());
    
    ASSERT_EQ(kcSq.size(), numModes);

    // Check with a tolerance (e.g., 1%)
    for (int i = 0; i < numModes; ++i) {
        ASSERT_NEAR(kcSq[i], kcSqAnalytical[i], kcSqAnalytical[i] * 0.01);
    }
}

// TODO: Add a test_parallel_eigenvalue_cutoff (Rule #25)
// This test would run the same problem with MPI and verify the results are identical.
```

**3. Implementation (TO PASS THE TEST):**
Create `src/hpcfem/physics/physics_waveguide_eigen.hpp` and `.cpp` (Rule \#17, \#16).

```cpp
// in src/hpcfem/physics/physics_waveguide_eigen.hpp
#pragma once
#include "mfem.hpp"
#include "hpcfem/core/physics_interface.hpp"
#include <vector>

namespace hpcfem {

/**
 * @brief Solves the 2D eigenvalue problem for waveguide modes.
 * @details This class solves the vector wave equation for the
 * transverse H-field (TE modes) or E-field (TM modes) to find
 * the cutoff wavenumbers (eigenvalues) and mode shapes (eigenvectors).
 * This is the foundation for port definitions.
 */
class PhysicsWaveguideEigen : public PhysicsInterface {
public:
    explicit PhysicsWaveguideEigen(mfem::Mesh& mesh);
    virtual ~PhysicsWaveguideEigen() = default;

    /**
     * @brief Assembles the FEM matrices for the eigenvalue problem.
     * @details Implements the weak form: (1/er * curl(H_t), curl(v)) = k_c^2 * (ur * H_t, v)
     */
    void setup() override;

    /**
     * @brief Solves the generalized eigenvalue problem [A]{x} = lambda [B]{x}.
     * @param numModes The number of eigenvalues to solve for.
     * @return A vector of eigenvalues (k_c^2).
     */
    std::vector<double> solve(int numModes);

    // TODO: Add methods to get the eigenvectors (mode shapes)

protected:
    void assembleA() override;
    void assembleB() override; // 'B' here is the mass matrix [M]

    // ND (Nedelec) collection for vector fields (Rule #20)
    mfem::ND_FECollection feCollection;
    mfem::FiniteElementSpace feSpace;

    // A = Stiffness Matrix (CurlCurl)
    mfem::BilinearForm a;
    // B = Mass Matrix
    mfem::BilinearForm b; 
    
    // TODO: Add support for SLEPc or ARPACK for eigenvalue solution
};

} // namespace hpcfem
```

  * **Implementation Note:** The `solve` method will interface with an eigensolver (like SLEPc, as used in MFEM's examples). You will assemble the `CurlCurl` integrator for matrix `a` and the `VectorFEMass` integrator for matrix `b`. The key is using `mfem::ND_FECollection`.

-----

## Step 2: The Driven Problem - Terminated Waveguide (ABC)

**1. Planning:**
Now we move to a 3D, *driven*, *complex-valued* problem. This is the vector Helmholtz equation: $\nabla \times (\frac{1}{\mu_r} \nabla \times \mathbf{E}) - k_0^2 \epsilon_r \mathbf{E} = 0$.

We will simulate a simple, straight rectangular waveguide.

  * **Source:** We'll use a simple current sheet to launch a wave.
  * **Walls:** PEC boundary conditions ($\mathbf{E} \times \mathbf{n} = 0$).
  * **Termination:** We will implement a simple Absorbing Boundary Condition (ABC) at the output port. This directly addresses the "terminated boundary" and "complex-value solving" problems.

**2. Test (WRITE THIS FIRST):**
Create a new file: `tests/test_physics_waveguide_driven.cpp`

```cpp
// in tests/test_physics_waveguide_driven.cpp
#include "gtest/gtest.h"
#include "mfem.hpp"
#include "hpcfem/physics/physics_waveguide_driven.hpp"
#include "hpcfem/core/fem_problem.hpp"

// Test for complex-valued driven problem with ABC termination
TEST(PhysicsWaveguideDriven, test_waveguide_abc_termination) {
    // 1. Setup the problem
    constexpr int dim = 3;
    constexpr int order = 1;
    std::string meshFile = "../testdata/testmesh_long_bar.mesh"; // (e.g., a 0.5x0.5x5.0 bar)
    
    constexpr double freq = 10e9; // 10 GHz
    constexpr double c0 = 299792458.0;
    constexpr double k0 = 2.0 * 3.1415926535 * freq / c0;

    hpcfem::FemProblem problem(dim, order, meshFile);
    
    // 2. Create and set the physics
    hpcfem::PhysicsWaveguideDriven physics(problem.getMesh(), k0);
    
    // 3. Define boundaries
    // Attribute 1: PEC Walls
    physics.addPecBoundary(1); 
    // Attribute 2: Source Plane (e.g., at z=0)
    physics.addSourceBoundary(2, 1.0); // 1.0 is dummy amplitude for now
    // Attribute 3: ABC Port (e.g., at z=5.0)
    physics.addAbcBoundary(3);

    problem.setPhysics(&physics);

    // 4. Solve the driven, complex-valued problem
    problem.solve(); // This will solve for the complex E-field

    // 5. Assertions
    mfem::GridFunction eField = problem.getSolution();
    
    // TODO: This test needs to be more robust.
    // A simple test: check the magnitude of the field at the end.
    // A better test:
    // 1. Get the E-field norm near the source.
    // 2. Get the E-field norm near the termination (ABC).
    // 3. Assert that norm_termination is significant (wave propagated).
    // 4. A proper test would compute S11 and assert it is low (e.g., < -20dB).
    // For now, we'll just check that the solve worked.
    ASSERT_GT(eField.Norml2(), 0.0);
}
```

**3. Implementation (TO PASS THE TEST):**
Create `src/hpcfem/physics/physics_waveguide_driven.hpp` and `.cpp`.

```cpp
// in src/hpcfem/physics/physics_waveguide_driven.hpp
#pragma once
#include "mfem.hpp"
#include "hpcfem/core/physics_interface.hpp"
#include <map>

namespace hpcfem {

/**
 * @brief Solves the 3D frequency-domain vector wave equation.
 * @details This class assembles the complex-valued system for:
 * Curl(1/ur * Curl(E)) - k0^2 * er * E = -j*k0*Z0*J
 * It handles PEC walls, current sources, and absorbing boundaries.
 */
class PhysicsWaveguideDriven : public PhysicsInterface {
public:
    PhysicsWaveguideDriven(mfem::Mesh& mesh, double k0);
    virtual ~PhysicsWaveguideDriven() = default;

    void setup() override;

    /** @brief Sets a PEC boundary condition (E x n = 0) */
    void addPecBoundary(int attribute);

    /** @brief Sets a first-order Absorbing Boundary Condition (ABC) */
    void addAbcBoundary(int attribute);

    /** @brief Sets a current source (J) */
    void addSourceBoundary(int attribute, double amplitude);

    // Overrides for complex-valued problem
    // We need to manage Real and Imaginary parts (Rule #1)
    
    // TODO: This interface needs to be updated to support complex
    // systems. A good "plain design" (Rule #2) is to have
    // the physics object own two solvers, one for Real and one for Imag.
    // Or, more simply, assemble the complex 2x2 block system:
    // [ K  -wM ] [Er] = [fr]
    // [ wM   K ] [Ei] = [fi]
    // Where K is CurlCurl and M is Mass.
    // This is a major implementation step.
    
    // TODO: For now, we will assume mfem::ComplexOperator
    // is available and we solve a complex system.
    
    // This requires mfem::ParComplexGridFunction, etc.
    mfem::GridFunction& getSolutionReal();
    mfem::GridFunction& getSolutionImag();


protected:
    void assembleA() override; // Assembles the full complex operator
    void assembleB() override; // Assembles the full complex RHS

    double k0; // Wavenumber
    mfem::ND_FECollection feCollection;
    mfem::FiniteElementSpace feSpace;

    // Boundary attributes
    mfem::Array<int> pecAttributes;
    mfem::Array<int> abcAttributes;
    mfem::Array<int> srcAttributes;
    
    // TODO: These will need to be complex-valued
    // mfem::BilinearForm a;
    // mfem::LinearForm b;
    // mfem::GridFunction eField;
};

} // namespace hpcfem
```

  * **Implementation Note:** This is the *hardest* step. You must implement the complex-valued system. The vector Helmholtz equation $\nabla \times (\nabla \times \mathbf{E}) - k_0^2 \mathbf{E} = 0$ is assembled as:
      * **Stiffness `[K]`:** `mfem::CurlCurlIntegrator`
      * **Mass `[M]`:** `mfem::VectorFEMassIntegrator`
      * **ABC `[C]`:** `mfem::VectorFEMassIntegrator` (on the boundary)
      * The full complex system is `[A] = [K] - k_0^2[M] + j k_0 [C]`. You will have to manage the real and imaginary parts of this system.

-----

## Step 3: The Advanced Port - Eigenfunction Expansion

**1. Planning:**
The ABC is an *approximation*. The *exact* solution is the hybrid FEM/Eigenfunction Expansion method (Chapter 10). This involves solving the 2D eigenproblem from Step 1, storing the modes, and using them as a boundary condition.

**2. Test (WRITE THIS FIRST):**
Create a new file: `tests/test_waveguide_port.cpp`

```cpp
// in tests/test_waveguide_port.cpp
#include "gtest/gtest.h"
#include "hpcfem/ports/waveguide_port.hpp"

TEST(WaveguidePort, test_port_mode_calculation) {
    // 1. Setup: Use the same 2D mesh as the eigenvalue test
    std::string meshFile = "../testdata/testmesh_rectangle_2d.mesh";
    constexpr double a = 2.286e-2;
    constexpr double b = 1.016e-2;
    constexpr int numModes = 5;

    // 2. Create the port object
    hpcfem::WaveguidePort port(meshFile, numModes, order);
    
    // 3. Action: Ask it to solve for its modes
    port.solveModes();

    // 4. Assertions
    // Check cutoff frequencies (same as Step 1 test)
    std::vector<double> kc = port.getCutoffWavenumbers();
    ASSERT_EQ(kc.size(), numModes);
    
    constexpr double PI = 3.141592653589793;
    double kc_TE10 = PI/a;
    ASSERT_NEAR(kc[0], kc_TE10, kc_TE10 * 0.01);

    // Check that it stored the eigenvectors (the mode shapes)
    // mfem::GridFunction* modeShape = port.getModeShape(0); // TE10
    // ASSERT_NE(modeShape, nullptr);
    // ASSERT_GT(modeShape->Norml2(), 0.0);
}

// Now, update the driven physics test
// in tests/test_physics_waveguide_driven.cpp
TEST(PhysicsWaveguideDriven, test_waveguide_eigenfunction_port) {
    // 1. Setup
    // ... same as before, but use a mesh of a discontinuity
    std::string meshFile = "../testdata/testmesh_waveguide_iris.mesh"; // TODO: Create this mesh
    
    // ...
    hpcfem::PhysicsWaveguideDriven physics(problem.getMesh(), k0);
    
    // 2. Create ports
    std::string portMeshFile = "../testdata/testmesh_rectangle_2d.mesh";
    hpcfem::WaveguidePort port1(portMeshFile, 5, order);
    hpcfem::WaveguidePort port2(portMeshFile, 5, order);

    // 3. Define boundaries
    physics.addPecBoundary(1); 
    // Attribute 2: Port 1 (Source)
    physics.addEigenfunctionPort(2, &port1, true); // true = isSource
    // Attribute 3: Port 2 (Termination)
    physics.addEigenfunctionPort(3, &port2, false); // false = isSource

    problem.setPhysics(&physics);
    problem.solve();

    // 4. Post-process for S-Parameters
    // TODO: Create SParameterCalculator class
    // hpcfem::SParameterCalculator sCalc;
    // double s11 = sCalc.getS11(physics, &port1);
    // double s21 = sCalc.getS21(physics, &port2);

    // 5. Assertions
    // For a simple through-waveguide (no iris), S11 ~ 0, S21 ~ 1
    // For an iris, we just check power conservation (Rule #10)
    // double powerConservation = std::pow(s11, 2) + std::pow(s21, 2);
    // ASSERT_NEAR(powerConservation, 1.0, 1e-3);
}
```

**3. Implementation (TO PASS THE TEST):**

**A. `src/hpcfem/ports/waveguide_port.hpp` (and `.cpp`)**

  * This class will hold a 2D `mfem::Mesh` of the port cross-section.
  * It will contain the `PhysicsWaveguideEigen` (from Step 1) as a private member.
  * `solveModes()` will call the eigensolver and store the `mfem::GridFunction` eigenvectors and `k_c` eigenvalues.
  * `getModeShape(int i)` will return the $i$-th mode.

**B. `src/hpcfem/physics/physics_waveguide_driven.hpp` (Updates)**

  * `addEigenfunctionPort(int attribute, WaveguidePort* port, bool isSource)`: This method is the core.
      * It will create a new custom `mfem::BilinearFormIntegrator` (e.g., `WaveguidePortIntegrator`).
      * This integrator will compute the port admittance matrix (using the modes from the `WaveguidePort` object) and add it to the global complex matrix `[A]`. This *is* the terminated boundary.
      * If `isSource == true`, it will use the `TE10` mode shape from the `port` object to build the *RHS vector* $\{b\}$, which excites the system.

**C. `src/hpcfem/analysis/s_parameter_calculator.hpp` (and `.cpp`)**

  * This is a post-processing class.
  * It will have methods like `getS11(PhysicsWaveguideDriven& physics, WaveguidePort& port)`.
  * This method will perform a complex-valued surface integral over the port's boundary attribute, projecting the total 3D E-field solution onto the 2D port modes.
  * `S_ij = \int_{Port_j} (\mathbf{E}_{total} \cdot \mathbf{e}_j) dA`
  * This calculation gives the complex amplitudes of the reflected ($S_{11}, S_{22}$) and transmitted ($S_{21}, S_{12}$) modes.

-----

## Step 4: Refactor, Document, and Stress Test

**1. Refactoring (Rule \#10, \#14):**

  * The `PhysicsWaveguideDriven` class will become large. Review it to ensure file-length (Rule \#14) and nesting (Rule \#4) rules are followed.
  * The complex-number handling (e.g., the 2x2 block matrix) should be clean and readable.

**2. Documenting (Rule \#5, \#6):**

  * Add all new classes (`PhysicsWaveguideEigen`, `PhysicsWaveguideDriven`, `WaveguidePort`, `SParameterCalculator`, `WaveguidePortIntegrator`) to `docs/hpcfem-doc/naming_registry.md`.
  * Ensure all public methods have Doxygen comments explaining their function.

**3. Stress Testing (Rule \#10):**

  * Modify the `test_waveguide_eigenfunction_port` test.
  * Load a much finer mesh (e.g., 1e6+ elements).
  * Run the simulation (in Release mode, `cmake-build-release`, Rule \#21) with `mpirun -n 16`.
  * The test passes if:
    1.  The parallel run completes without crashing (Rule \#25).
    2.  The S-parameter results are consistent with the low-resolution serial run.
    3.  The solver (e.g., Hypre AMG, `SolverHypreAmg`) (Rule \#12) converges efficiently.