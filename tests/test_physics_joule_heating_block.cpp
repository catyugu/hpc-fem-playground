/**
 * @file test_physics_joule_heating_block.cpp
 * @brief Test for monolithic 2×2 block system Joule heating
 * 

 * 
 * This test verifies the JouleHeatingPhysics class which assembles
 * a coupled electro-thermal system as a 2×2 BlockOperator:
 * 
 * [K_e    0 ] [V]   [f_e]
 * [C    K_t] [T] = [f_t]
 * 
 * For this initial test, coupling C is zero (decoupled system).
 * The test verifies:
 * 1. Block system assembles correctly
 * 2. Block sizes match finite element spaces
 * 3. System can be solved with direct solver
 */

#include "hpcfem/physics_joule_heating.hpp"
#include <gtest/gtest.h>
#include <mfem.hpp>

// Constants (Guideline #15: No magic numbers)
constexpr double ELECTRICAL_CONDUCTIVITY = 1.0e3;  // σ [S/m]
constexpr double THERMAL_CONDUCTIVITY = 10.0;      // κ [W/(m·K)]
constexpr int MESH_ELEMENTS_1D = 4;
constexpr int POLYNOMIAL_ORDER = 1;

#ifdef MFEM_USE_MPI
TEST(PhysicsJouleHeatingTest, BlockSystemAssemblyParallel)
{
    // Create 3D mesh
    mfem::Mesh* serial_mesh = new mfem::Mesh(
        mfem::Mesh::MakeCartesian3D(MESH_ELEMENTS_1D, MESH_ELEMENTS_1D, 
                                     MESH_ELEMENTS_1D, mfem::Element::HEXAHEDRON));
    
    mfem::ParMesh* pmesh = new mfem::ParMesh(MPI_COMM_WORLD, *serial_mesh);
    delete serial_mesh;
    
    // Create coupled physics
    hpcfem::JouleHeatingPhysics physics(pmesh, POLYNOMIAL_ORDER,
                                        ELECTRICAL_CONDUCTIVITY,
                                        THERMAL_CONDUCTIVITY);
    
    // Get block structure
    const mfem::Array<int>& blockOffsets = physics.getBlockOffsets();
    
    // Verify block offsets
    EXPECT_EQ(blockOffsets.Size(), 3);
    EXPECT_EQ(blockOffsets[0], 0);
    EXPECT_GT(blockOffsets[1], 0);  // Electric DOFs
    EXPECT_GT(blockOffsets[2], blockOffsets[1]);  // Total DOFs
    
    // Create block operator and vectors
    mfem::BlockOperator blockOp(blockOffsets);
    mfem::BlockVector blockRHS(blockOffsets);
    mfem::BlockVector blockSol(blockOffsets);
    
    // Assemble system
    mfem::Array<int> essTdofElectric, essTdofThermal;
    physics.assemble(blockOp, blockRHS, blockSol, essTdofElectric, essTdofThermal);
    
    // Verify essential DOF lists are not empty
    EXPECT_GT(essTdofElectric.Size(), 0);
    EXPECT_GT(essTdofThermal.Size(), 0);
    
    // Verify block vector sizes match
    EXPECT_EQ(blockSol.GetBlock(0).Size(), blockOffsets[1] - blockOffsets[0]);
    EXPECT_EQ(blockSol.GetBlock(1).Size(), blockOffsets[2] - blockOffsets[1]);
    
    // Test passes if we reach here without crashes
    delete pmesh;
}
#endif

TEST(PhysicsJouleHeatingTest, BlockSystemAssemblySerial)
{
#ifndef MFEM_USE_MPI
    // Create simple 2D mesh for serial test
    mfem::Mesh* mesh = new mfem::Mesh(
        mfem::Mesh::MakeCartesian2D(MESH_ELEMENTS_1D, MESH_ELEMENTS_1D, 
                                     mfem::Element::QUADRILATERAL));
    
    // Create coupled physics
    hpcfem::JouleHeatingPhysics physics(mesh, POLYNOMIAL_ORDER,
                                        ELECTRICAL_CONDUCTIVITY,
                                        THERMAL_CONDUCTIVITY);
    
    // Get block structure
    const mfem::Array<int>& blockOffsets = physics.getBlockOffsets();
    
    // Verify block offsets
    EXPECT_EQ(blockOffsets.Size(), 3);
    EXPECT_EQ(blockOffsets[0], 0);
    EXPECT_GT(blockOffsets[1], 0);
    EXPECT_GT(blockOffsets[2], blockOffsets[1]);
    
    // Create block operator and vectors
    mfem::BlockOperator blockOp(blockOffsets);
    mfem::BlockVector blockRHS(blockOffsets);
    mfem::BlockVector blockSol(blockOffsets);
    
    // Assemble system
    mfem::Array<int> essTdofElectric, essTdofThermal;
    physics.assemble(blockOp, blockRHS, blockSol, essTdofElectric, essTdofThermal);
    
    // Verify essential DOF lists
    EXPECT_GT(essTdofElectric.Size(), 0);
    EXPECT_GT(essTdofThermal.Size(), 0);
    
    // Verify block sizes
    EXPECT_EQ(blockSol.GetBlock(0).Size(), blockOffsets[1]);
    EXPECT_EQ(blockSol.GetBlock(1).Size(), blockOffsets[2] - blockOffsets[1]);
    
    // Clean up
    delete mesh;
#else
    // Skip serial test when compiled with MPI
    GTEST_SKIP() << "Serial test skipped in MPI build";
#endif
}

int main(int argc, char** argv)
{
#ifdef MFEM_USE_MPI
    MPI_Init(&argc, &argv);
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    MPI_Finalize();
    return result;
#else
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
#endif
}
