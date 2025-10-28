/**
 * @file solver_hypre_amg.hpp
 * @brief Concrete HYPRE BoomerAMG solver implementation
 * 
 * Phase: 2, Step: 2.2
 * 
 * This class wraps MFEM's HypreBoomerAMG preconditioner as a complete
 * linear solver that implements the SolverInterface.
 */

#ifndef HPCFEM_SOLVER_HYPRE_AMG_HPP
#define HPCFEM_SOLVER_HYPRE_AMG_HPP

#include "solver_interface.hpp"
#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class HypreAmgSolver
 * @brief HYPRE BoomerAMG-based linear solver
 * 
 * Encapsulates the setup and usage of HYPRE's BoomerAMG algebraic
 * multigrid preconditioner with a conjugate gradient (CG) solver.
 */
class HypreAmgSolver : public SolverInterface
{
public:
    /**
     * @brief Constructor with solver parameters
     * 
     * @param relTol Relative tolerance for convergence (default: 1e-12)
     * @param maxIter Maximum number of iterations (default: 2000)
     * @param printLevel Verbosity level (default: 0 = silent)
     */
    HypreAmgSolver(double relTol = 1.0e-12, 
                   int maxIter = 2000, 
                   int printLevel = 0);

    /**
     * @brief Destructor - cleans up solver resources
     */
    ~HypreAmgSolver() override = default;

    /**
     * @brief Solve the linear system Ax = b using AMG-preconditioned CG
     * 
     * @param A System matrix (HypreParMatrix for parallel)
     * @param b Right-hand side vector
     * @param x Solution vector (input: initial guess, output: solution)
     */
#ifdef MFEM_USE_MPI
    void solve(const mfem::HypreParMatrix& A,
              const mfem::Vector& b,
              mfem::Vector& x) override;
#else
    void solve(const mfem::SparseMatrix& A,
              const mfem::Vector& b,
              mfem::Vector& x) override;
#endif

    /**
     * @brief Get the number of iterations from the last solve
     * @return Number of iterations
     */
    int getNumIterations() const { return numIterations_; }

    /**
     * @brief Get the final residual norm from the last solve
     * @return Final residual norm
     */
    double getFinalNorm() const { return finalNorm_; }

private:
    double relTol_;        ///< Relative convergence tolerance
    int maxIter_;          ///< Maximum iterations
    int printLevel_;       ///< Verbosity level
    int numIterations_;    ///< Actual iterations from last solve
    double finalNorm_;     ///< Final residual norm from last solve
};

} // namespace hpcfem

#endif // HPCFEM_SOLVER_HYPRE_AMG_HPP
