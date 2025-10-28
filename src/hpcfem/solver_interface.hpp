/**
 * @file solver_interface.hpp
 * @brief Abstract base class for all linear solvers in the hpcfem library
 * 

 * 
 * This interface defines the contract for all linear solvers.
 * Concrete implementations will wrap various solver backends (HYPRE, PETSc, etc.)
 */

#ifndef HPCFEM_SOLVER_INTERFACE_HPP
#define HPCFEM_SOLVER_INTERFACE_HPP

#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class SolverInterface
 * @brief Pure abstract base class for linear solvers
 * 
 * All solver implementations must inherit from this class and implement
 * the solve() method. This enables polymorphic solver selection at runtime.
 */
class SolverInterface
{
public:
    /**
     * @brief Virtual destructor for proper cleanup
     */
    virtual ~SolverInterface() = default;

    /**
     * @brief Solve the linear system Ax = b
     * 
     * @param A System matrix (must be assembled and finalized)
     * @param b Right-hand side vector
     * @param x Solution vector (input: initial guess, output: solution)
     * 
     * @note This is a pure virtual function - all derived classes must implement it
     */
#ifdef MFEM_USE_MPI
    virtual void solve(const mfem::HypreParMatrix& A, 
                      const mfem::Vector& b, 
                      mfem::Vector& x) = 0;
#else
    virtual void solve(const mfem::SparseMatrix& A,
                      const mfem::Vector& b,
                      mfem::Vector& x) = 0;
#endif
};

} // namespace hpcfem

#endif // HPCFEM_SOLVER_INTERFACE_HPP
