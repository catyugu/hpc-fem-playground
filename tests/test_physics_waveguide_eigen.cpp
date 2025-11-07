/**
 * @file test_physics_cavity_eigen.cpp
 * @brief Test for 3D Maxwell eigenvalue problem solver
 * 
 * This test validates the 3D electromagnetic eigenvalue solver by computing
 * cavity modes. It uses a simple rectangular cavity geometry where analytical
 * solutions are known.
 * 
 * Follows TDD principles (Rule #10) and tests both serial and parallel (Rule #25)
 */

#include "gtest/gtest.h"
#include "mfem.hpp"
#include "hpcfem/physics/physics_cavity_eigen.hpp"

#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

// Constants (Rule #15: No magic numbers)
constexpr double PI = 3.141592653589793;

// WR-90 standard rectangular waveguide dimensions
constexpr double CAVITY_LENGTH_X = 0.02286;  // 22.86 mm (0.9 inches) - width 'a'
constexpr double CAVITY_LENGTH_Y = 0.01016;  // 10.16 mm (0.4 inches) - height 'b'
constexpr double CAVITY_LENGTH_Z = 0.05000;  // 50 mm - cavity length (arbitrary for testing)

constexpr int NUM_MODES_TO_COMPUTE = 10;
constexpr int POLYNOMIAL_ORDER = 1;
constexpr double EIGENVALUE_TOLERANCE = 0.10; // 10% tolerance (3D FEM is less accurate than analytical)

/**
 * @brief Test 3D rectangular cavity eigenvalue calculation
 * @details Validates the eigenvalue solver against analytical solutions
 * for a WR-90 rectangular cavity resonator. For a cavity with PEC walls,
 * the eigenvalues are: λ = k² = (mπ/Lx)² + (nπ/Ly)² + (pπ/Lz)²
 * where m,n,p are mode indices (integers ≥ 0, but not all zero)
 * 
 * WR-90 dimensions: 22.86mm × 10.16mm (standard waveguide)
 */
TEST(PhysicsCavityEigen, test_3d_cavity_modes)
{
    // 1. Setup the problem - create a 3D rectangular cavity mesh
    constexpr int elementsX = 20;
    constexpr int elementsY = 12;
    constexpr int elementsZ = 40;
    
#ifdef MFEM_USE_MPI
    // Create 3D Cartesian mesh
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian3D(
        elementsX, elementsY, elementsZ, mfem::Element::HEXAHEDRON,
        CAVITY_LENGTH_X, CAVITY_LENGTH_Y, CAVITY_LENGTH_Z);
    mfem::ParMesh mesh(MPI_COMM_WORLD, serialMesh);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        std::cout << "\n=== 3D Cavity Eigenvalue Test (WR-90) ===" << std::endl;
        std::cout << "Cavity dimensions: " << CAVITY_LENGTH_X*1000 << " x " 
                  << CAVITY_LENGTH_Y*1000 << " x " << CAVITY_LENGTH_Z*1000 << " mm" << std::endl;
        std::cout << "Mesh: " << elementsX << " x " << elementsY << " x " 
                  << elementsZ << " elements" << std::endl;
    }
#else
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
        elementsX, elementsY, elementsZ, mfem::Element::HEXAHEDRON,
        CAVITY_LENGTH_X, CAVITY_LENGTH_Y, CAVITY_LENGTH_Z);
        
    std::cout << "\n=== 3D Cavity Eigenvalue Test (WR-90) ===" << std::endl;
    std::cout << "Cavity dimensions: " << CAVITY_LENGTH_X*1000 << " x " 
              << CAVITY_LENGTH_Y*1000 << " x " << CAVITY_LENGTH_Z*1000 << " mm" << std::endl;
    std::cout << "Mesh: " << elementsX << " x " << elementsY << " x " 
              << elementsZ << " elements" << std::endl;
#endif

    // 2. Create the physics object
    hpcfem::PhysicsCavityEigen physics(&mesh, POLYNOMIAL_ORDER);
    
    // 3. Solve the eigenvalue problem
    std::vector<double> eigenvalues = physics.solveEigenvalues(NUM_MODES_TO_COMPUTE);

    // 4. Compute analytical Maxwell cavity solutions for lowest modes
    // For a rectangular PEC cavity (dimensions Lx,Ly,Lz) modal eigenvalues
    // are k^2 = (m*pi/Lx)^2 + (n*pi/Ly)^2 + (p*pi/Lz)^2 where (m,n,p)
    // are integers. TM modes require m,n,p >= 1. TE modes allow zeros but
    // not the trivial (0,0,0). We'll generate both TE and TM sets up to a
    // modest index, deduplicate, sort, and match FEM eigenvalues to the
    // nearest analytical values.
    std::vector<double> analyticalEigenvalues;
    const int maxIndex = 4; // increase if you need more analytical modes
    const double pi = PI;

    for (int m = 0; m <= maxIndex; ++m)
    {
        for (int n = 0; n <= maxIndex; ++n)
        {
            for (int p = 0; p <= maxIndex; ++p)
            {
                if (m == 0 && n == 0 && p == 0) { continue; }

                // TE: any non-trivial combination
                double k2_te = std::pow(m * pi / CAVITY_LENGTH_X, 2) +
                               std::pow(n * pi / CAVITY_LENGTH_Y, 2) +
                               std::pow(p * pi / CAVITY_LENGTH_Z, 2);
                analyticalEigenvalues.push_back(k2_te);

                // TM: require all indices positive
                if (m >= 1 && n >= 1 && p >= 1)
                {
                    double k2_tm = k2_te; // same formula, different mode type
                    analyticalEigenvalues.push_back(k2_tm);
                }
            }
        }
    }

    // Deduplicate and sort analytical values
    std::sort(analyticalEigenvalues.begin(), analyticalEigenvalues.end());
    const double dedup_eps = 1e-8;
    std::vector<double> uniqueAnalytical;
    for (double v : analyticalEigenvalues)
    {
        if (uniqueAnalytical.empty() || std::abs(v - uniqueAnalytical.back()) > dedup_eps)
        {
            uniqueAnalytical.push_back(v);
        }
    }
    analyticalEigenvalues.swap(uniqueAnalytical);

    ASSERT_EQ(eigenvalues.size(), NUM_MODES_TO_COMPUTE)
        << "Expected " << NUM_MODES_TO_COMPUTE << " eigenvalues";

    // Match FEM eigenvalues to nearest analytical eigenvalues (one-to-one)
    std::vector<bool> used(analyticalEigenvalues.size(), false);
    for (size_t i = 0; i < eigenvalues.size(); ++i)
    {
        double fem = eigenvalues[i];
        size_t best_j = 0;
        double best_diff = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j < analyticalEigenvalues.size(); ++j)
        {
            if (used[j]) { continue; }
            double diff = std::abs(fem - analyticalEigenvalues[j]);
            if (diff < best_diff)
            {
                best_diff = diff;
                best_j = j;
            }
        }
        double rel_err = best_diff / analyticalEigenvalues[best_j];

#ifdef MFEM_USE_MPI
        if (rank == 0)
        {
            std::cout << "Match FEM mode " << i << ": FEM=" << fem
                      << " analytical=" << analyticalEigenvalues[best_j]
                      << " rel_err=" << rel_err*100.0 << "%\n";
        }
#else
        std::cout << "Match FEM mode " << i << ": FEM=" << fem
                  << " analytical=" << analyticalEigenvalues[best_j]
                  << " rel_err=" << rel_err*100.0 << "%\n";
#endif

        EXPECT_LE(rel_err, EIGENVALUE_TOLERANCE) << "Mode " << i << " mismatch: rel_err=" << rel_err;
        used[best_j] = true;
    }
}

/**
 * @brief Test that eigenvectors are properly computed
 * @details Validates that mode shapes can be extracted as grid functions
 */
TEST(PhysicsCavityEigen, test_mode_shapes)
{
    // Setup
    constexpr int elementsX = 3;
    constexpr int elementsY = 2;
    constexpr int elementsZ = 2;
    
#ifdef MFEM_USE_MPI
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian3D(
        elementsX, elementsY, elementsZ, mfem::Element::HEXAHEDRON,
        0.05, 0.03, 0.04);
    mfem::ParMesh mesh(MPI_COMM_WORLD, serialMesh);
#else
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(
        elementsX, elementsY, elementsZ, mfem::Element::HEXAHEDRON,
        0.05, 0.03, 0.04);
#endif

    hpcfem::PhysicsCavityEigen physics(&mesh, POLYNOMIAL_ORDER);
    
    // Solve
    constexpr int numModesTest = 3;
    std::vector<double> eigenvalues = physics.solveEigenvalues(numModesTest);
    
    // Get eigenvectors
    std::vector<mfem::Vector> eigenvectors = physics.getEigenvectors();
    
    ASSERT_EQ(eigenvectors.size(), numModesTest) 
        << "Expected " << numModesTest << " eigenvectors";
    
    // Check that eigenvectors are non-zero
    for (size_t i = 0; i < eigenvectors.size(); ++i)
    {
        double norm = eigenvectors[i].Norml2();
        EXPECT_GT(norm, 0.0) << "Eigenvector " << i << " should be non-zero";
    }
    
    // Check that we can get mode shapes as grid functions
    for (int i = 0; i < numModesTest; ++i)
    {
#ifdef MFEM_USE_MPI
        mfem::ParGridFunction* modeShape = physics.getModeShape(i);
#else
        mfem::GridFunction* modeShape = physics.getModeShape(i);
#endif
        ASSERT_NE(modeShape, nullptr) << "Mode shape " << i << " should not be null";
        EXPECT_GT(modeShape->Norml2(), 0.0) << "Mode shape " << i << " should be non-zero";
        delete modeShape;
    }
}

// MPI initialization for parallel tests (Rule #25)
int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    
#ifdef MFEM_USE_MPI
    MPI_Init(&argc, &argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
#else
    return RUN_ALL_TESTS();
#endif
}
