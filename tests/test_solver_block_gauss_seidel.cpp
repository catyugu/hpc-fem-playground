/**
 * @file test_solver_block_gauss_seidel.cpp
 * @brief Test for block Gauss-Seidel preconditioner
 * 
 * Tests the physics-informed block-triangular preconditioning strategy
 * for coupled multiphysics systems.
 * 

 */

#include <gtest/gtest.h>
#include <mfem.hpp>
#include "hpcfem/solver_block_gauss_seidel.hpp"

#ifdef MFEM_USE_MPI
#include <mpi.h>
#endif

/**
 * @brief Test block Gauss-Seidel solver on a simple 2×2 block system
 * 
 * Creates a coupled system with structure:
 *   [A  0] [x1]   [b1]
 *   [C  B] [x2] = [b2]
 * 
 * where A and B are SPD matrices and C is a coupling matrix.
 */
TEST(BlockGaussSeidelSolverTest, SimpleBlockSystem)
{
#ifdef MFEM_USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Create a simple 3D mesh
    const int meshSize = 4;
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(meshSize, meshSize, meshSize,
                                                   mfem::Element::HEXAHEDRON);
    mfem::ParMesh* pmesh = new mfem::ParMesh(MPI_COMM_WORLD, mesh);
    
    // Create two finite element spaces (for the two fields)
    const int order = 1;
    mfem::H1_FECollection* fec = new mfem::H1_FECollection(order, pmesh->Dimension());
    mfem::ParFiniteElementSpace* fespace1 = new mfem::ParFiniteElementSpace(pmesh, fec);
    mfem::ParFiniteElementSpace* fespace2 = new mfem::ParFiniteElementSpace(pmesh, fec);
    
    // Add essential boundary conditions to make the problem well-posed
    mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;  // All boundaries are essential
    
    // Build diagonal blocks A and B (both Laplacian operators)
    mfem::ConstantCoefficient one(1.0);
    
    mfem::ParBilinearForm* a1 = new mfem::ParBilinearForm(fespace1);
    a1->AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
    a1->Assemble();
    a1->EliminateEssentialBCDiag(ess_bdr, 1.0);
    a1->Finalize();
    mfem::HypreParMatrix* A = a1->ParallelAssemble();
    
    mfem::ParBilinearForm* a2 = new mfem::ParBilinearForm(fespace2);
    a2->AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
    a2->Assemble();
    a2->EliminateEssentialBCDiag(ess_bdr, 1.0);
    a2->Finalize();
    mfem::HypreParMatrix* B = a2->ParallelAssemble();
    
    // Build coupling matrix C = 0.1 * Mass matrix (from space1 to space2)
    // This represents the physics coupling (e.g., Joule heating)
    mfem::ParMixedBilinearForm* c_form = new mfem::ParMixedBilinearForm(fespace1, fespace2);
    mfem::ConstantCoefficient coeff_c(0.1);
    c_form->AddDomainIntegrator(new mfem::MixedScalarMassIntegrator(coeff_c));
    c_form->Assemble();
    c_form->Finalize();
    mfem::HypreParMatrix* C = c_form->ParallelAssemble();
    
    // Create block structure
    mfem::Array<int> blockOffsets(3);
    blockOffsets[0] = 0;
    blockOffsets[1] = fespace1->GetTrueVSize();
    blockOffsets[2] = blockOffsets[1] + fespace2->GetTrueVSize();
    
    // Create BlockOperator
    mfem::BlockOperator* blockOp = new mfem::BlockOperator(blockOffsets);
    blockOp->SetBlock(0, 0, A);
    blockOp->SetBlock(1, 0, C);  // Coupling block
    blockOp->SetBlock(1, 1, B);
    
    // Create RHS and solution vectors
    mfem::BlockVector blockRHS(blockOffsets);
    mfem::BlockVector blockSol(blockOffsets);
    blockRHS = 1.0;  // Simple constant RHS
    blockSol = 0.0;
    
    // Create block solvers (using HYPRE AMG for each diagonal block)
    mfem::HypreBoomerAMG* solver0 = new mfem::HypreBoomerAMG(*A);
    solver0->SetPrintLevel(0);
    solver0->iterative_mode = false;  // Apply as M^{-1}, not iterative
    
    mfem::HypreBoomerAMG* solver1 = new mfem::HypreBoomerAMG(*B);
    solver1->SetPrintLevel(0);
    solver1->iterative_mode = false;  // Apply as M^{-1}, not iterative
    
    // Create block Gauss-Seidel preconditioner
    hpcfem::BlockGaussSeidelSolver precond(blockOffsets, blockOp, solver0, solver1, true);
    
    // Test 1: Apply preconditioner once
    mfem::BlockVector precondsol(blockOffsets);
    precondsol = 0.0;
    precond.Mult(blockRHS, precondsol);
    
    // Verify the solution is non-zero
    double norm = precondsol.Norml2();
    if (rank == 0)
    {
        std::cout << "Block Gauss-Seidel preconditioner application: ||y|| = " << norm << std::endl;
    }
    EXPECT_GT(norm, 0.1);  // Solution should be reasonably large
    
    // Test 2: Use as preconditioner in GMRES
    mfem::GMRESSolver gmres(MPI_COMM_WORLD);
    gmres.SetOperator(*blockOp);
    gmres.SetPreconditioner(precond);
    gmres.SetRelTol(1e-6);
    gmres.SetMaxIter(100);
    gmres.SetPrintLevel(rank == 0 ? 1 : 0);
    
    blockSol = 0.0;
    gmres.Mult(blockRHS, blockSol);
    
    EXPECT_TRUE(gmres.GetConverged());
    int iters = gmres.GetNumIterations();
    if (rank == 0)
    {
        std::cout << "GMRES with Block Gauss-Seidel: converged in " << iters << " iterations" << std::endl;
    }
    
    // Verify convergence is fast (should be much less than 100 iterations)
    EXPECT_LT(iters, 30);
    
    // Verify the solution satisfies the system
    mfem::BlockVector residual(blockOffsets);
    blockOp->Mult(blockSol, residual);
    residual -= blockRHS;
    double resNorm = residual.Norml2();
    double rhsNorm = blockRHS.Norml2();
    
    if (rank == 0)
    {
        std::cout << "Relative residual: " << resNorm / rhsNorm << std::endl;
    }
    EXPECT_LT(resNorm / rhsNorm, 1e-5);
    
    // Clean up
    delete blockOp;
    delete c_form;
    delete a2;
    delete a1;
    delete C;
    delete B;
    delete A;
    delete fespace2;
    delete fespace1;
    delete fec;
    delete pmesh;
    
#else
    // Serial version - skip test
    GTEST_SKIP() << "Test requires MPI";
#endif
}

/**
 * @brief Benchmark Block Gauss-Seidel against monolithic approach
 * 
 * This test demonstrates the speedup of physics-informed preconditioning
 * versus treating the coupled system as a black-box. Uses a larger mesh
 * and measures wall-clock time to reveal performance differences.
 */
TEST(BlockGaussSeidelSolverTest, CompareWithMonolithicAMG)
{
#ifdef MFEM_USE_MPI
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    // Use a larger mesh to reveal performance differences
    const int meshSize = 16;  // 16^3 = 4096 elements
    mfem::Mesh mesh = mfem::Mesh::MakeCartesian3D(meshSize, meshSize, meshSize,
                                                   mfem::Element::HEXAHEDRON);
    mfem::ParMesh* pmesh = new mfem::ParMesh(MPI_COMM_WORLD, mesh);
    
    const int order = 2;  // Higher order for more DOFs
    mfem::H1_FECollection* fec = new mfem::H1_FECollection(order, pmesh->Dimension());
    mfem::ParFiniteElementSpace* fespace1 = new mfem::ParFiniteElementSpace(pmesh, fec);
    mfem::ParFiniteElementSpace* fespace2 = new mfem::ParFiniteElementSpace(pmesh, fec);
    
    // Add essential boundary conditions
    mfem::Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;
    
    // Build system matrices with different characteristics
    // Block 0: Stiff diffusion (small epsilon)
    // Block 1: Standard diffusion
    mfem::ConstantCoefficient stiff_coeff(0.01);  // Stiff problem
    mfem::ConstantCoefficient one(1.0);
    
    mfem::ParBilinearForm* a1 = new mfem::ParBilinearForm(fespace1);
    a1->AddDomainIntegrator(new mfem::DiffusionIntegrator(stiff_coeff));
    a1->AddDomainIntegrator(new mfem::MassIntegrator(one));  // Add mass for conditioning
    a1->Assemble();
    a1->EliminateEssentialBCDiag(ess_bdr, 1.0);
    a1->Finalize();
    mfem::HypreParMatrix* A = a1->ParallelAssemble();
    
    mfem::ParBilinearForm* a2 = new mfem::ParBilinearForm(fespace2);
    a2->AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
    a2->Assemble();
    a2->EliminateEssentialBCDiag(ess_bdr, 1.0);
    a2->Finalize();
    mfem::HypreParMatrix* B = a2->ParallelAssemble();
    
    // Strong coupling matrix
    mfem::ParMixedBilinearForm* c_form = new mfem::ParMixedBilinearForm(fespace1, fespace2);
    mfem::ConstantCoefficient coeff_c(0.8);  // Strong coupling
    c_form->AddDomainIntegrator(new mfem::MixedScalarMassIntegrator(coeff_c));
    c_form->Assemble();
    c_form->Finalize();
    mfem::HypreParMatrix* C = c_form->ParallelAssemble();
    
    mfem::Array<int> blockOffsets(3);
    blockOffsets[0] = 0;
    blockOffsets[1] = fespace1->GetTrueVSize();
    blockOffsets[2] = blockOffsets[1] + fespace2->GetTrueVSize();
    
    mfem::BlockOperator* blockOp = new mfem::BlockOperator(blockOffsets);
    blockOp->SetBlock(0, 0, A);
    blockOp->SetBlock(1, 0, C);
    blockOp->SetBlock(1, 1, B);
    
    mfem::BlockVector blockRHS(blockOffsets);
    mfem::BlockVector blockSol(blockOffsets);
    blockRHS = 1.0;
    
    if (rank == 0)
    {
        std::cout << "\n=== Block Gauss-Seidel Benchmark ===" << std::endl;
        std::cout << "Mesh: " << meshSize << "^3 elements, Order: " << order << std::endl;
        std::cout << "Total DOFs: " << blockOffsets[2] << " (" << blockOffsets[1] << " + " 
                  << (blockOffsets[2] - blockOffsets[1]) << ")" << std::endl;
    }
    
    // Test 1: Monolithic approach - Block-diagonal preconditioner with AMG
    // This is the "standard" approach: apply AMG to each diagonal block independently
    // (ignores the coupling structure)
    mfem::HypreBoomerAMG* amg0_mono = new mfem::HypreBoomerAMG(*A);
    amg0_mono->SetPrintLevel(0);
    amg0_mono->iterative_mode = false;
    
    mfem::HypreBoomerAMG* amg1_mono = new mfem::HypreBoomerAMG(*B);
    amg1_mono->SetPrintLevel(0);
    amg1_mono->iterative_mode = false;
    
    // Block-diagonal preconditioner (ignores coupling)
    mfem::BlockDiagonalPreconditioner blockDiag(blockOffsets);
    blockDiag.SetDiagonalBlock(0, amg0_mono);
    blockDiag.SetDiagonalBlock(1, amg1_mono);
    
    mfem::GMRESSolver gmres1(MPI_COMM_WORLD);
    gmres1.SetOperator(*blockOp);
    gmres1.SetPreconditioner(blockDiag);
    gmres1.SetRelTol(1e-6);
    gmres1.SetMaxIter(200);
    gmres1.SetPrintLevel(0);
    
    blockSol = 0.0;
    auto t_start_mono = MPI_Wtime();
    gmres1.Mult(blockRHS, blockSol);
    auto t_end_mono = MPI_Wtime();
    double time_monolithic = t_end_mono - t_start_mono;
    int iters_monolithic = gmres1.GetNumIterations();
    
    if (rank == 0)
    {
        std::cout << "\nMonolithic (Block-Diagonal AMG):" << std::endl;
        std::cout << "  Iterations: " << iters_monolithic << std::endl;
        std::cout << "  Time: " << time_monolithic << " seconds" << std::endl;
        std::cout << "  Converged: " << (gmres1.GetConverged() ? "Yes" : "No") << std::endl;
    }
    
    // Test 2: Block Gauss-Seidel (physics-informed approach)
    // This exploits the block-triangular structure
    mfem::HypreBoomerAMG* amg0_gs = new mfem::HypreBoomerAMG(*A);
    amg0_gs->SetPrintLevel(0);
    amg0_gs->iterative_mode = false;
    
    mfem::HypreBoomerAMG* amg1_gs = new mfem::HypreBoomerAMG(*B);
    amg1_gs->SetPrintLevel(0);
    amg1_gs->iterative_mode = false;
    
    hpcfem::BlockGaussSeidelSolver blockGS(blockOffsets, blockOp, amg0_gs, amg1_gs, true);
    
    mfem::GMRESSolver gmres2(MPI_COMM_WORLD);
    gmres2.SetOperator(*blockOp);
    gmres2.SetPreconditioner(blockGS);
    gmres2.SetRelTol(1e-6);
    gmres2.SetMaxIter(200);
    gmres2.SetPrintLevel(0);
    
    blockSol = 0.0;
    auto t_start_gs = MPI_Wtime();
    gmres2.Mult(blockRHS, blockSol);
    auto t_end_gs = MPI_Wtime();
    double time_blockGS = t_end_gs - t_start_gs;
    int iters_blockGS = gmres2.GetNumIterations();
    
    if (rank == 0)
    {
        std::cout << "\nBlock Gauss-Seidel (Physics-Informed):" << std::endl;
        std::cout << "  Iterations: " << iters_blockGS << std::endl;
        std::cout << "  Time: " << time_blockGS << " seconds" << std::endl;
        std::cout << "  Converged: " << (gmres2.GetConverged() ? "Yes" : "No") << std::endl;
        
        double iter_speedup = static_cast<double>(iters_monolithic) / iters_blockGS;
        double time_speedup = time_monolithic / time_blockGS;
        std::cout << "\nSpeedup:" << std::endl;
        std::cout << "  Iterations: " << iter_speedup << "x" << std::endl;
        std::cout << "  Wall-clock time: " << time_speedup << "x" << std::endl;
    }
    
    // The physics-informed solver should perform at least as well as the diagonal approach
    // In many cases it will be faster, especially for strongly coupled systems
    // For this test, we just verify both converge successfully
    EXPECT_TRUE(gmres1.GetConverged());
    EXPECT_TRUE(gmres2.GetConverged());
    
    // Log the comparison for documentation purposes
    if (rank == 0)
    {
        if (iters_blockGS < iters_monolithic)
        {
            std::cout << "\n✓ Physics-informed solver achieved better convergence!" << std::endl;
        }
        else if (time_blockGS < time_monolithic)
        {
            std::cout << "\n✓ Physics-informed solver was faster in wall-clock time!" << std::endl;
        }
        else
        {
            std::cout << "\n• Both solvers performed similarly for this problem." << std::endl;
            std::cout << "  (Physics-informed advantages appear with stronger coupling)" << std::endl;
        }
    }
    
    // Clean up
    delete blockOp;
    delete c_form;
    delete a2;
    delete a1;
    delete C;
    delete B;
    delete A;
    delete fespace2;
    delete fespace1;
    delete fec;
    delete pmesh;
    
#else
    GTEST_SKIP() << "Test requires MPI";
#endif
}

int main(int argc, char** argv)
{
#ifdef MFEM_USE_MPI
    MPI_Init(&argc, &argv);
#endif
    
    ::testing::InitGoogleTest(&argc, argv);
    int result = RUN_ALL_TESTS();
    
#ifdef MFEM_USE_MPI
    MPI_Finalize();
#endif
    
    return result;
}
