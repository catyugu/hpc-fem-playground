/**
 * @file solver_block_gauss_seidel.hpp
 * @brief Block Gauss-Seidel preconditioner for coupled multiphysics systems
 * 
 * Implements physics-informed block-triangular preconditioning strategy that
 * exploits the one-way coupling structure in systems like Joule heating.
 * 

 */

#ifndef HPCFEM_SOLVER_BLOCK_GAUSS_SEIDEL_HPP
#define HPCFEM_SOLVER_BLOCK_GAUSS_SEIDEL_HPP

#include <mfem.hpp>

namespace hpcfem
{

/**
 * @class BlockGaussSeidelSolver
 * @brief Physics-informed block-triangular preconditioner
 * 
 * For a 2×2 block system with one-way coupling:
 *   [K_e  0 ] [V]   [f_e]
 *   [C   K_t] [T] = [f_t]
 * 
 * The block Gauss-Seidel iteration performs:
 *   1. Solve K_e * V = f_e           (electrical problem)
 *   2. Update RHS: f_t' = f_t - C*V  (coupling term)
 *   3. Solve K_t * T = f_t'          (thermal problem)
 * 
 * This is "physics-informed" because it respects the physical coupling
 * structure (electricity drives heat, but not vice versa).
 * 
 * Expected performance: 4× faster iteration count vs monolithic AMG.
 */
class BlockGaussSeidelSolver : public mfem::Solver
{
public:
    /**
     * @brief Constructor
     * @param blockOffsets Array defining block structure [0, n1, n1+n2]
     * @param blockOperator The 2×2 block system operator
     * @param solver0 Solver for diagonal block (0,0) - typically K_e
     * @param solver1 Solver for diagonal block (1,1) - typically K_t
     * @param ownSolvers If true, this class owns and will delete the solvers
     * 
     * Note: blockOperator must have block (1,0) set (the coupling matrix C).
     *       Block (0,1) should be nullptr (one-way coupling).
     */
    BlockGaussSeidelSolver(const mfem::Array<int>& blockOffsets,
                          mfem::BlockOperator* blockOperator,
                          mfem::Solver* solver0,
                          mfem::Solver* solver1,
                          bool ownSolvers = false);
    
    /**
     * @brief Destructor
     */
    ~BlockGaussSeidelSolver() override;
    
    /**
     * @brief Apply the block Gauss-Seidel preconditioner
     * @param x Input vector (RHS), output vector (solution increment)
     * 
     * Performs one sweep of block forward substitution.
     */
    void Mult(const mfem::Vector& x, mfem::Vector& y) const override;
    
    /**
     * @brief Set the operator (required by mfem::Solver interface)
     * @param op Operator to set (should be a BlockOperator)
     * 
     * Note: This updates the internal reference but does not take ownership.
     */
    void SetOperator(const mfem::Operator& op) override;

private:
    // Block structure
    mfem::Array<int> blockOffsets_;  ///< Block offsets [0, n1, n1+n2]
    
    // Block operator (not owned)
    mfem::BlockOperator* blockOp_;
    
    // Diagonal block solvers (may be owned)
    mfem::Solver* solver0_;  ///< Solver for block (0,0)
    mfem::Solver* solver1_;  ///< Solver for block (1,1)
    bool ownSolvers_;
    
    // Coupling matrix (not owned, reference from blockOp_)
    mfem::Operator* couplingOp_;  ///< Block (1,0) coupling operator C
    
    // Work vectors (mutable for use in const Mult)
    mutable mfem::BlockVector xBlock_;  ///< Block view of input
    mutable mfem::BlockVector yBlock_;  ///< Block view of output
    mutable mfem::Vector tmpVec_;       ///< Temporary for coupling term
};

} // namespace hpcfem

#endif // HPCFEM_SOLVER_BLOCK_GAUSS_SEIDEL_HPP
