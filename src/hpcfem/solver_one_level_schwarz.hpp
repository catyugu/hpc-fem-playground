/**
 * @file solver_one_level_schwarz.hpp
 * @brief One-level overlapping Schwarz (Additive Schwarz) preconditioner
 * 
 * Phase: 2, Step: 2.1
 * 
 * This class implements the classical one-level additive Schwarz method:
 * 
 *   M^{-1} = ∑_{i=1}^{N} R_i^T A_i^{-1} R_i
 * 
 * where:
 *   - N is the number of subdomains (MPI processes)
 *   - R_i is the restriction operator to subdomain i
 *   - A_i is the local matrix on subdomain i
 * 
 * This method is known to have poor scalability - iteration count increases
 * with the number of subdomains. This will be demonstrated in benchmarks.
 */

#ifndef HPCFEM_SOLVER_ONE_LEVEL_SCHWARZ_HPP
#define HPCFEM_SOLVER_ONE_LEVEL_SCHWARZ_HPP

#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class OneLevelSchwarz
 * @brief One-level overlapping Schwarz preconditioner
 * 
 * Implements the additive Schwarz preconditioner as an mfem::Solver.
 * This preconditioner performs local solves on each subdomain (MPI process)
 * and combines results using restriction/prolongation operators.
 * 
 * The local solver uses HYPRE BoomerAMG for efficiency.
 */
class OneLevelSchwarz : public mfem::Solver
{
public:
    /**
     * @brief Constructor - sets up the preconditioner
     * 
     * @param A Global parallel matrix (HypreParMatrix)
     * 
     * This constructor extracts the local subdomain matrix and sets up
     * the local solver (AMG) for each MPI process.
     */
    explicit OneLevelSchwarz(const mfem::HypreParMatrix& A);

    /**
     * @brief Destructor - cleans up local solver resources
     */
    ~OneLevelSchwarz() override;

    /**
     * @brief Apply the preconditioner: y = M^{-1} x
     * 
     * @param x Input vector
     * @param y Output vector (preconditioned result)
     * 
     * This performs the additive Schwarz operation:
     * 1. Restrict x to local subdomain
     * 2. Solve local system A_i y_i = x_i
     * 3. Prolongate y_i back to global vector
     * 4. Sum contributions from all subdomains
     */
    void Mult(const mfem::Vector& x, mfem::Vector& y) const override;

    /**
     * @brief Set the operator (required by mfem::Solver interface)
     * 
     * @param op Operator to set
     * 
     * Note: For this preconditioner, the operator is set in the constructor.
     * This method is provided for interface compliance but does nothing.
     */
    void SetOperator(const mfem::Operator& op) override;

private:
    const mfem::HypreParMatrix* globalMatrix_;  ///< Reference to global matrix
    mfem::SparseMatrix* localMatrix_;           ///< Local subdomain matrix (diagonal block)
    mfem::Solver* localSolver_;                 ///< Local solver (GSSmoother)
    
    mutable mfem::Vector localX_;  ///< Temporary vector for local operations
    mutable mfem::Vector localY_;  ///< Temporary vector for local operations
    
    int localSize_;        ///< Size of local subdomain problem
};

} // namespace hpcfem

#endif // HPCFEM_SOLVER_ONE_LEVEL_SCHWARZ_HPP
