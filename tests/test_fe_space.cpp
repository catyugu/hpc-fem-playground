#include <gtest/gtest.h>
#include "mfem.hpp"
#include "hpcfem/core/finite_element_space.hpp"
 
#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

TEST(FiniteElementSpaceTests, H1SpaceCreationSerial)
{
    std::string mesh_path = "testdata/testmesh_cube.mesh";
    mfem::Mesh mesh(mesh_path, 1, 1);
    const int order = 1;
    const int vdim = 1;

#ifdef MFEM_USE_MPI
    // Ensure MPI is initialized when constructing ParMesh
    int _mpi_initialized = 0;
    MPI_Initialized(&_mpi_initialized);
    bool _mpi_init_by_test = false;
    if (!_mpi_initialized) { MPI_Init(NULL, NULL); _mpi_init_by_test = true; }

    // Construct ParMesh from serial mesh for MPI builds
    mfem::ParMesh *pmesh = new mfem::ParMesh(MPI_COMM_WORLD, mesh);
    hpcfem::FiniteElementSpace fes(pmesh, order, vdim, hpcfem::FiniteElementSpace::Type::H1);
    auto* mfemSpace = fes.getMfemSpace();
    ASSERT_NE(mfemSpace, nullptr);
    EXPECT_EQ(mfemSpace->GetVDim(), vdim);
    delete pmesh;
    if (_mpi_init_by_test) { MPI_Finalize(); }
#else
    hpcfem::FiniteElementSpace fes(&mesh, order, vdim, hpcfem::FiniteElementSpace::Type::H1);
    auto* mfemSpace = fes.getMfemSpace();

    ASSERT_NE(mfemSpace, nullptr);
    EXPECT_EQ(mfemSpace->GetVDim(), vdim);
#endif
}
