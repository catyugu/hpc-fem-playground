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
 * @brief Schwarz domain decomposition solver using PCG acceleration
 * 
 * This solver combines domain decomposition with Krylov acceleration:
 * - Uses Additive Schwarz preconditioner (local AMG solves)
 * - Accelerated by Preconditioned Conjugate Gradient (PCG)
 * - Each subdomain solved independently with HYPRE BoomerAMG
 * - CG provides optimal convergence for SPD systems
 * 
 * Key advantages over Richardson iteration:
 * - Superlinear convergence (vs linear for Richardson)
 * - No parameter tuning required (no omega)
 * - Typically 2-3x fewer iterations
 * 
 * This is particularly effective for:
 * - Large-scale parallel problems
 * - Problems with good domain decomposition structure
 * - SPD (symmetric positive definite) systems
 */
class DdmSchwarzSolver : public SolverInterface
{
public:
    /**
     * @brief Constructor for DDM Schwarz solver with PCG
     * @param relTol Relative tolerance for PCG convergence
     * @param maxIter Maximum number of CG iterations
     * @param printLevel Print level (0=none, 1=summary, 2+=detailed)
     * @param omega Unused (kept for interface compatibility)
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
     * @brief Solve the parallel linear system using DDM-PCG
     * @param A System matrix (parallel)
     * @param b Right-hand side vector
     * @param x Solution vector (input: initial guess, output: solution)
     * 
     * Implementation uses MFEM's CGSolver with HypreBoomerAMG preconditioner.
     * The AMG hierarchy is built once and reused for all CG iterations.
     */
    void solve(const mfem::HypreParMatrix& A,
              const mfem::Vector& b,
              mfem::Vector& x) override;
#else
    /**
     * @brief Solve the serial linear system using CG with GS preconditioner
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
