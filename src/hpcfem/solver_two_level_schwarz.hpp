/**
 * @file solver_two_level_schwarz.hpp
 * @brief Two-level overlapping Schwarz preconditioner with AMG coarse solver
 * 
 * Phase: 2, Step: 2.2
 * 
 * This class implements the two-level additive Schwarz method with coarse-grid correction:
 * 
 *   M^{-1} = Ψ A_H^{-1} Ψ^T + ∑_{i=1}^{N} R_i^T A_i^{-1} R_i
 * 
 * where:
 *   - First term: coarse-grid correction using Galerkin projection A_H = Ψ^T A Ψ
 *   - Second term: local subdomain solves (same as one-level method)
 *   - Ψ is the prolongation operator from coarse to fine grid
 *   - A_H^{-1} is solved using HYPRE BoomerAMG
 * 
 * This method provides scalable performance - iteration count remains nearly constant
 * as the number of subdomains increases, unlike the one-level method.
 */

#ifndef HPCFEM_SOLVER_TWO_LEVEL_SCHWARZ_HPP
#define HPCFEM_SOLVER_TWO_LEVEL_SCHWARZ_HPP

#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class TwoLevelSchwarz
 * @brief Two-level overlapping Schwarz preconditioner with AMG coarse solver
 * 
 * Implements the two-level Schwarz preconditioner combining local subdomain
 * solves with a global coarse-grid correction. This provides scalability
 * by enabling fast global information propagation through the coarse grid.
 */
class TwoLevelSchwarz : public mfem::Solver
{
public:
    /**
     * @brief Constructor - sets up the two-level preconditioner
     * 
     * @param A Global parallel matrix (HypreParMatrix)
     * @param fespace Pointer to the fine-grid finite element space
     * 
     * This constructor:
     * 1. Extracts local subdomain matrices for local solves
     * 2. Constructs coarse finite element space
     * 3. Builds coarse-grid matrix using Galerkin projection
     * 4. Sets up AMG solver for coarse-grid problem
     */
    TwoLevelSchwarz(const mfem::HypreParMatrix& A, 
                    mfem::ParFiniteElementSpace* fespace);

    /**
     * @brief Destructor - cleans up solver resources
     */
    ~TwoLevelSchwarz() override;

    /**
     * @brief Apply the preconditioner: y = M^{-1} x
     * 
     * @param x Input vector
     * @param y Output vector (preconditioned result)
     * 
     * This performs the two-level Schwarz operation:
     * 1. Apply local subdomain solves: y_local = ∑ R_i^T A_i^{-1} R_i x
     * 2. Apply coarse-grid correction: y_coarse = Ψ A_H^{-1} Ψ^T x
     * 3. Combine: y = y_local + y_coarse
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
    // Fine-grid components (for local solves)
    const mfem::HypreParMatrix* globalMatrix_;  ///< Reference to fine-grid matrix
    mfem::SparseMatrix* localMatrix_;           ///< Local subdomain matrix
    mfem::Solver* localSolver_;                 ///< Local GS solver
    
    // Coarse-grid components
    mfem::ParFiniteElementSpace* fineSpace_;    ///< Fine-grid FE space (not owned)
    mfem::ParFiniteElementSpace* coarseSpace_;  ///< Coarse-grid FE space (owned)
    mfem::H1_FECollection* coarseFec_;          ///< Coarse-grid FE collection (owned)
    
    mfem::HypreParMatrix* coarseMatrix_;        ///< Coarse-grid matrix A_H (owned)
    mfem::HypreSolver* coarseSolver_;           ///< AMG solver for A_H (owned)
    
    mfem::HypreParMatrix* prolongation_;        ///< Prolongation Ψ: coarse → fine (owned)
    
    // Temporary vectors for local operations
    mutable mfem::Vector localX_;   ///< Local input vector
    mutable mfem::Vector localY_;   ///< Local output vector
    
    // Temporary vectors for coarse operations
    mutable mfem::Vector coarseX_;  ///< Coarse input vector
    mutable mfem::Vector coarseY_;  ///< Coarse output vector
    mutable mfem::Vector tempFine_; ///< Temporary fine vector for restriction
    
    int localSize_;   ///< Size of local subdomain problem
    int coarseSize_;  ///< Size of coarse-grid problem
};

} // namespace hpcfem

#endif // HPCFEM_SOLVER_TWO_LEVEL_SCHWARZ_HPP
