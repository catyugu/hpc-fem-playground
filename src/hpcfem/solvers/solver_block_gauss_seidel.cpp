/**
 * @file solvers/solver_block_gauss_seidel.cpp
 * @brief Implementation of block Gauss-Seidel preconditioner
 */

#include "hpcfem/solvers/solver_block_gauss_seidel.hpp"

namespace hpcfem
{

BlockGaussSeidelSolver::BlockGaussSeidelSolver(const mfem::Array<int>& blockOffsets,
                                               mfem::BlockOperator* blockOperator,
                                               mfem::Solver* solver0,
                                               mfem::Solver* solver1,
                                               bool ownSolvers)
    : mfem::Solver(blockOffsets[blockOffsets.Size()-1]),
      blockOffsets_(blockOffsets),
      blockOp_(blockOperator),
      solver0_(solver0),
      solver1_(solver1),
      ownSolvers_(ownSolvers),
      couplingOp_(nullptr),
      xBlock_(blockOffsets),
      yBlock_(blockOffsets),
      tmpVec_(blockOffsets[2] - blockOffsets[1])
{
    couplingOp_ = &blockOp_->GetBlock(1, 0);
    if (!couplingOp_)
    {
        mfem::mfem_error("BlockGaussSeidelSolver: Block (1,0) coupling operator is required!");
    }
    solver0_->SetOperator(blockOp_->GetBlock(0, 0));
    solver1_->SetOperator(blockOp_->GetBlock(1, 1));
}

BlockGaussSeidelSolver::~BlockGaussSeidelSolver()
{
    if (ownSolvers_)
    {
        delete solver0_;
        delete solver1_;
    }
}

void BlockGaussSeidelSolver::Mult(const mfem::Vector& x, mfem::Vector& y) const
{
    xBlock_.Update(const_cast<double*>(x.GetData()), blockOffsets_);
    yBlock_.Update(y.GetData(), blockOffsets_);

    mfem::Vector x0, x1;
    xBlock_.GetBlockView(0, x0);
    xBlock_.GetBlockView(1, x1);

    mfem::Vector y0, y1;
    yBlock_.GetBlockView(0, y0);
    yBlock_.GetBlockView(1, y1);

    y0 = 0.0;
    solver0_->Mult(x0, y0);

    couplingOp_->Mult(y0, tmpVec_);

    y1 = 0.0;
    mfem::Vector rhs1(x1.Size());
    subtract(x1, tmpVec_, rhs1);
    solver1_->Mult(rhs1, y1);
}

void BlockGaussSeidelSolver::SetOperator(const mfem::Operator& op)
{
    blockOp_ = const_cast<mfem::BlockOperator*>(dynamic_cast<const mfem::BlockOperator*>(&op));
    if (!blockOp_)
    {
        mfem::mfem_error("BlockGaussSeidelSolver::SetOperator: Operator must be a BlockOperator!");
    }
    couplingOp_ = &blockOp_->GetBlock(1, 0);
    if (!couplingOp_)
    {
        mfem::mfem_error("BlockGaussSeidelSolver::SetOperator: Block (1,0) coupling operator is required!");
    }
    solver0_->SetOperator(blockOp_->GetBlock(0, 0));
    solver1_->SetOperator(blockOp_->GetBlock(1, 1));
}

} // namespace hpcfem
