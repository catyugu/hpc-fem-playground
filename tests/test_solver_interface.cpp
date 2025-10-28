/**
 * @file test_solver_interface.cpp
 * @brief Compile-only test for the SolverInterface abstract base class
 * 
 * Phase: 2, Step: 2.1
 * 
 * This test verifies that:
 * 1. The header can be included without errors
 * 2. A mock solver can inherit from SolverInterface
 * 3. The interface compiles correctly
 */

#include "hpcfem/solver_interface.hpp"
#include <iostream>
#include <cstdlib>

using namespace hpcfem;

/**
 * @class MockSolver
 * @brief A minimal concrete implementation for testing
 */
class MockSolver : public SolverInterface
{
public:
    MockSolver() = default;
    ~MockSolver() override = default;

#ifdef MFEM_USE_MPI
    void solve(const mfem::HypreParMatrix& A,
              const mfem::Vector& b,
              mfem::Vector& x) override
    {
        // Trivial implementation for testing
        x = b;
    }
#else
    void solve(const mfem::SparseMatrix& A,
              const mfem::Vector& b,
              mfem::Vector& x) override
    {
        // Trivial implementation for testing
        x = b;
    }
#endif
};

int main(int argc, char *argv[])
{
    // Initialize MPI if available
#ifdef MFEM_USE_MPI
    mfem::MPI_Session mpi(argc, argv);
    int myid = mpi.WorldRank();
#else
    int myid = 0;
#endif

    if (myid == 0)
    {
        std::cout << "Testing SolverInterface..." << std::endl;
    }

    // Test 1: Can create a mock solver
    MockSolver mock_solver;

    // Test 2: Can use polymorphic pointer
    SolverInterface* solver_ptr = &mock_solver;

    // Test 3: Basic type checking
    if (solver_ptr != nullptr)
    {
        if (myid == 0)
        {
            std::cout << "Test passed: SolverInterface compiles and supports polymorphism" << std::endl;
        }
        return EXIT_SUCCESS;
    }

    return EXIT_FAILURE;
}
