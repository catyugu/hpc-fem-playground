/**
 * @file test_solver_interface_mock.cpp
 * @brief Mock test for SolverInterface contract verification
 * 
 * Phase: 1, Step: 1.2
 * 
 * This test verifies the contract of the SolverInterface:
 * - Mock solver can be instantiated
 * - solve() method can be called
 * - Interface polymorphism works correctly
 */

#include "hpcfem/solver_interface.hpp"
#include <gtest/gtest.h>

using namespace hpcfem;

/**
 * @class MockSolver
 * @brief Mock implementation of SolverInterface for testing
 */
class MockSolver : public SolverInterface
{
public:
    MockSolver() : solveCallCount(0) {}
    ~MockSolver() override = default;

#ifdef MFEM_USE_MPI
    void solve(const mfem::HypreParMatrix& A,
              const mfem::Vector& b,
              mfem::Vector& x) override
    {
        solveCallCount++;
        // Mock implementation: just copy b to x
        x = b;
    }
#else
    void solve(const mfem::SparseMatrix& A,
              const mfem::Vector& b,
              mfem::Vector& x) override
    {
        solveCallCount++;
        // Mock implementation: just copy b to x
        x = b;
    }
#endif

    int getSolveCallCount() const { return solveCallCount; }

private:
    int solveCallCount;
};

/**
 * @brief Test that MockSolver implements SolverInterface correctly
 */
TEST(SolverInterfaceTest, MockSolverInstantiation)
{
    MockSolver solver;
    EXPECT_EQ(solver.getSolveCallCount(), 0);
}

/**
 * @brief Test that solve() method can be called through interface
 */
TEST(SolverInterfaceTest, SolveMethodCall)
{
    MockSolver solver;
    
    // Create a simple 2x2 matrix and vectors
    constexpr int SIZE = 2;
    
#ifdef MFEM_USE_MPI
    // For parallel build, we need a simple test matrix
    // Create a minimal problem to get a valid HypreParMatrix
    MPI_Comm comm = MPI_COMM_WORLD;
    
    // Create a simple 1D mesh
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian1D(SIZE);
    mfem::ParMesh pmesh(comm, serialMesh);
    
    // Create H1 space
    mfem::H1_FECollection fec(1, pmesh.Dimension());
    mfem::ParFiniteElementSpace fespace(&pmesh, &fec);
    
    // Create bilinear form to get a matrix
    mfem::ParBilinearForm a(&fespace);
    mfem::ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
    a.Assemble();
    a.Finalize();
    
    mfem::HypreParMatrix* Aptr = a.ParallelAssemble();
    mfem::HypreParMatrix& A = *Aptr;
#else
    // For serial build, create a simple sparse matrix
    mfem::SparseMatrix A(SIZE, SIZE);
    A.Set(0, 0, 1.0);
    A.Set(1, 1, 1.0);
    A.Finalize();
#endif
    
    mfem::Vector b(SIZE);
    b = 1.0;
    
    mfem::Vector x(SIZE);
    x = 0.0;
    
    // Call solve through interface
    SolverInterface* solverPtr = &solver;
    solverPtr->solve(A, b, x);
    
    // Verify solve was called
    EXPECT_EQ(solver.getSolveCallCount(), 1);
    
    // Verify mock behavior (x should equal b)
    EXPECT_DOUBLE_EQ(x(0), b(0));
    EXPECT_DOUBLE_EQ(x(1), b(1));
    
#ifdef MFEM_USE_MPI
    delete Aptr;
#endif
}

/**
 * @brief Test polymorphism through interface pointer
 */
TEST(SolverInterfaceTest, Polymorphism)
{
    std::unique_ptr<SolverInterface> solver = std::make_unique<MockSolver>();
    
    EXPECT_NE(solver.get(), nullptr);
    
    // Verify we can call solve through base class pointer
    constexpr int SIZE = 2;
    
#ifdef MFEM_USE_MPI
    MPI_Comm comm = MPI_COMM_WORLD;
    
    // Create a simple 1D mesh
    mfem::Mesh serialMesh = mfem::Mesh::MakeCartesian1D(SIZE);
    mfem::ParMesh pmesh(comm, serialMesh);
    
    // Create H1 space
    mfem::H1_FECollection fec(1, pmesh.Dimension());
    mfem::ParFiniteElementSpace fespace(&pmesh, &fec);
    
    // Create bilinear form to get a matrix
    mfem::ParBilinearForm a(&fespace);
    mfem::ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
    a.Assemble();
    a.Finalize();
    
    mfem::HypreParMatrix* Aptr = a.ParallelAssemble();
    mfem::HypreParMatrix& A = *Aptr;
#else
    mfem::SparseMatrix A(SIZE, SIZE);
    A.Set(0, 0, 1.0);
    A.Set(1, 1, 1.0);
    A.Finalize();
#endif
    
    mfem::Vector b(SIZE);
    b = 2.0;
    mfem::Vector x(SIZE);
    x = 0.0;
    
    solver->solve(A, b, x);
    
    EXPECT_DOUBLE_EQ(x(0), 2.0);
    EXPECT_DOUBLE_EQ(x(1), 2.0);
    
#ifdef MFEM_USE_MPI
    delete Aptr;
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
