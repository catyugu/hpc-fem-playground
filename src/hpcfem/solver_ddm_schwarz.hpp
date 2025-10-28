/**
 * @file solver_ddm_schwarz.hpp
 * @brief Domain Decomposition Method (DDM) with Schwarz preconditioner
 * 
 * This solver implements an overlapping Schwarz domain decomposition method
 * for parallel solution of linear systems. Each MPI process solves its local
 * subdomain problem, and the solutions are combined through iterative refinement.
 */

#ifndef HPCFEM_SOLVER_DDM_SCHWARZ_HPP
#define HPCFEM_SOLVER_DDM_SCHWARZ_HPP

#include "solver_interface.hpp"
#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class DdmSchwarzSolver
 * @brief Schwarz domain decomposition solver for parallel systems
 * 
 * This solver uses an additive Schwarz preconditioned Richardson iteration:
 * 1. Each subdomain solves its local problem independently (using AMG)
 * 2. Local solutions are combined with relaxation
 * 3. Ghost/halo values are exchanged via MPI
 * 4. Iterate until global convergence
 * 
 * This is particularly effective for:
 * - Large-scale parallel problems
 * - Problems with good domain decomposition structure
 * - When communication costs need to be minimized
 */
class DdmSchwarzSolver : public SolverInterface
{
public:
    /**
     * @brief Constructor for DDM Schwarz solver
     * @param relTol Relative tolerance for convergence
     * @param maxIter Maximum number of outer iterations
     * @param printLevel Print level (0=none, 1=summary, 2+=detailed)
     * @param omega Relaxation parameter (default 1.0 for standard Richardson)
     */
    DdmSchwarzSolver(double relTol = 1.0e-6,
                     int maxIter = 100,
                     int printLevel = 0,
                     double omega = 1.0);

    /**
     * @brief Destructor
     */
    ~DdmSchwarzSolver() override;

#ifdef MFEM_USE_MPI
    /**
     * @brief Solve the parallel linear system using DDM
     * @param A System matrix (parallel)
     * @param b Right-hand side vector
     * @param x Solution vector (input: initial guess, output: solution)
     */
    void solve(const mfem::HypreParMatrix& A,
              const mfem::Vector& b,
              mfem::Vector& x) override;
#else
    /**
     * @brief Solve the serial linear system (falls back to direct AMG)
     * @param A System matrix (serial)
     * @param b Right-hand side vector
     * @param x Solution vector (input: initial guess, output: solution)
     */
    void solve(const mfem::SparseMatrix& A,
              const mfem::Vector& b,
              mfem::Vector& x) override;
#endif

    /**
     * @brief Get the number of iterations performed
     * @return Iteration count
     */
    int getNumIterations() const;

    /**
     * @brief Get the final residual norm
     * @return Final norm
     */
    double getFinalNorm() const;

private:
    double relTol_;
    int maxIter_;
    int printLevel_;
    double omega_;
    
    int numIterations_;
    double finalNorm_;
};

} // namespace hpcfem

#endif // HPCFEM_SOLVER_DDM_SCHWARZ_HPP
