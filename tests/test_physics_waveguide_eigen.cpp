/**
 * @file test_physics_waveguide_eigen.cpp
 * @brief Unit tests for 2D waveguide eigenmode solver (TE/TM)
 */

#include "gtest/gtest.h"
#include "mfem.hpp"
#include "hpcfem/physics/physics_waveguide_eigen.hpp"

#include <vector>
#include <cmath>

constexpr double PI = 3.141592653589793;
constexpr double WR90_A = 0.02286; // width (a) in meters
constexpr double WR90_B = 0.01016; // height (b)

// Helper: analytical k_c^2 for rectangular waveguide
double analytical_kc2(int m, int n, double a, double b)
{
    double kx = (m * PI / a);
    double ky = (n * PI / b);
    return kx*kx + ky*ky;
}

TEST(PhysicsWaveguideEigen, TM_modes_rectangular)
{
    const int nx = 40, ny = 20;
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, 1, WR90_A, WR90_B);
#ifdef MFEM_USE_MPI
    mfem::ParMesh mesh(MPI_COMM_WORLD, serialMesh);
    hpcfem::PhysicsWaveguideEigen physics(&mesh, 1, hpcfem::WaveguideModeType::TM);
#else
    mfem::Mesh mesh = std::move(serialMesh);
    hpcfem::PhysicsWaveguideEigen physics(&mesh, 1, hpcfem::WaveguideModeType::TM);
#endif
    const int numModes = 6;
    std::vector<double> eigen = physics.solveEigenvalues(numModes);
    std::cout << "Computed eigenvalues:" << std::endl;
    for (size_t i=0;i<eigen.size();++i) std::cout << "  ["<<i<<"] "<<eigen[i]<<std::endl;

    ASSERT_EQ((int)eigen.size(), numModes);

    // Compare first few modes against analytical (m,n) ordering
    // lowest non-zero TM has m>=1,n>=1 so TM11 is first
    double target_tm11 = analytical_kc2(1,1,WR90_A,WR90_B);
    // find nearest computed
    double best = 1e300; int best_i = -1;
    for (int i=0;i<numModes;++i){ double d = std::abs(eigen[i]-target_tm11); if (d<best){best=d;best_i=i;} }
    double rel_err = std::abs(eigen[best_i]-target_tm11)/target_tm11;
    EXPECT_LE(rel_err, 0.10) << "TM11 relative error too large";
}

TEST(PhysicsWaveguideEigen, TE_modes_rectangular)
{
    const int nx = 40, ny = 20;
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian2D(nx, ny, mfem::Element::QUADRILATERAL, 1, WR90_A, WR90_B);
#ifdef MFEM_USE_MPI
    mfem::ParMesh mesh(MPI_COMM_WORLD, serialMesh);
    hpcfem::PhysicsWaveguideEigen physics(&mesh, 1, hpcfem::WaveguideModeType::TE);
#else
    mfem::Mesh mesh = std::move(serialMesh);
    hpcfem::PhysicsWaveguideEigen physics(&mesh, 1, hpcfem::WaveguideModeType::TE);
#endif
    const int numModes = 6;
    std::vector<double> eigen = physics.solveEigenvalues(numModes);

    ASSERT_EQ((int)eigen.size(), numModes);

    // TE10 is the lowest TE mode
    double target_te10 = analytical_kc2(1,0,WR90_A,WR90_B);
    double best = 1e300; int best_i = -1;
    for (int i=0;i<numModes;++i){ double d = std::abs(eigen[i]-target_te10); if (d<best){best=d;best_i=i;} }
    double rel_err = std::abs(eigen[best_i]-target_te10)/target_te10;
    EXPECT_LE(rel_err, 0.10) << "TE10 relative error too large";
}

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
