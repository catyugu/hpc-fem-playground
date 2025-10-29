/**
 * @file solvers/solver_block_gauss_seidel.hpp
 * @brief Block Gauss-Seidel preconditioner for coupled multiphysics systems
 */

#ifndef HPCFEM_SOLVER_BLOCK_GAUSS_SEIDEL_HPP
#define HPCFEM_SOLVER_BLOCK_GAUSS_SEIDEL_HPP

#include <mfem.hpp>

namespace hpcfem
{

class BlockGaussSeidelSolver : public mfem::Solver
{
public:
    BlockGaussSeidelSolver(const mfem::Array<int>& blockOffsets,
                          mfem::BlockOperator* blockOperator,
                          mfem::Solver* solver0,
                          mfem::Solver* solver1,
                          bool ownSolvers = false);

    ~BlockGaussSeidelSolver() override;

    void Mult(const mfem::Vector& x, mfem::Vector& y) const override;

    void SetOperator(const mfem::Operator& op) override;

private:
    mfem::Array<int> blockOffsets_;
    mfem::BlockOperator* blockOp_;
    mfem::Solver* solver0_;
    mfem::Solver* solver1_;
    bool ownSolvers_;
    mfem::Operator* couplingOp_;
    mutable mfem::BlockVector xBlock_;
    mutable mfem::BlockVector yBlock_;
    mutable mfem::Vector tmpVec_;
};

} // namespace hpcfem

#endif // HPCFEM_SOLVER_BLOCK_GAUSS_SEIDEL_HPP
