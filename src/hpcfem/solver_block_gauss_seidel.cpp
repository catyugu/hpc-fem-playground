/**
 * @file solver_block_gauss_seidel.cpp
 * @brief Implementation of block Gauss-Seidel preconditioner
 * 
 * Phase: 3, Step: 3.2
 */

#include "solver_block_gauss_seidel.hpp"

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
      tmpVec_(blockOffsets[2] - blockOffsets[1])  // Size of block 1
{
    // Get the coupling operator from block (1,0)
    couplingOp_ = &blockOp_->GetBlock(1, 0);
    
    // Verify it's not null (we need the coupling for Gauss-Seidel)
    if (!couplingOp_)
    {
        mfem::mfem_error("BlockGaussSeidelSolver: Block (1,0) coupling operator is required!");
    }
    
    // Set the operators for the block solvers
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
    // Create block views for input and output vectors
    xBlock_.Update(const_cast<double*>(x.GetData()), blockOffsets_);
    yBlock_.Update(y.GetData(), blockOffsets_);
    
    // Extract block vectors
    mfem::Vector x0, x1;
    xBlock_.GetBlockView(0, x0);
    xBlock_.GetBlockView(1, x1);
    
    mfem::Vector y0, y1;
    yBlock_.GetBlockView(0, y0);
    yBlock_.GetBlockView(1, y1);
    
    // Step 1: Solve K_e * y0 = x0 (electrical problem)
    y0 = 0.0;
    solver0_->Mult(x0, y0);
    
    // Step 2: Compute coupling term: tmpVec = C * y0
    couplingOp_->Mult(y0, tmpVec_);
    
    // Step 3: Solve K_t * y1 = (x1 - tmpVec) (thermal problem with coupling)
    y1 = 0.0;
    mfem::Vector rhs1(x1.Size());
    subtract(x1, tmpVec_, rhs1);
    solver1_->Mult(rhs1, y1);
}

void BlockGaussSeidelSolver::SetOperator(const mfem::Operator& op)
{
    // Cast to BlockOperator
    blockOp_ = const_cast<mfem::BlockOperator*>(dynamic_cast<const mfem::BlockOperator*>(&op));
    
    if (!blockOp_)
    {
        mfem::mfem_error("BlockGaussSeidelSolver::SetOperator: Operator must be a BlockOperator!");
    }
    
    // Update coupling operator reference
    couplingOp_ = &blockOp_->GetBlock(1, 0);
    
    if (!couplingOp_)
    {
        mfem::mfem_error("BlockGaussSeidelSolver::SetOperator: Block (1,0) coupling operator is required!");
    }
    
    // Update operators for block solvers
    solver0_->SetOperator(blockOp_->GetBlock(0, 0));
    solver1_->SetOperator(blockOp_->GetBlock(1, 1));
}

} // namespace hpcfem
