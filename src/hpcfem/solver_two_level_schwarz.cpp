/**
 * @file solver_two_level_schwarz.cpp
 * @brief Implementation of two-level overlapping Schwarz preconditioner
 * 
 * Phase: 2, Step: 2.2
 */

#include "hpcfem/solver_two_level_schwarz.hpp"

namespace hpcfem
{

TwoLevelSchwarz::TwoLevelSchwarz(const mfem::HypreParMatrix& A, 
                                 mfem::ParFiniteElementSpace* fespace)
    : mfem::Solver(A.Height()),
      globalMatrix_(&A),
      localMatrix_(nullptr),
      localSolver_(nullptr),
      fineSpace_(fespace),
      coarseSpace_(nullptr),
      coarseFec_(nullptr),
      coarseMatrix_(nullptr),
      coarseSolver_(nullptr),
      prolongation_(nullptr)
{
    // ========== Part 1: Local subdomain setup (same as one-level) ==========
    
    // Get local subdomain matrix (diagonal block)
    localMatrix_ = new mfem::SparseMatrix();
    A.GetDiag(*localMatrix_);
    localSize_ = localMatrix_->Height();
    
    // Allocate temporary vectors for local operations
    localX_.SetSize(localSize_);
    localY_.SetSize(localSize_);
    
    // Create local Gauss-Seidel smoother
    localSolver_ = new mfem::GSSmoother(*localMatrix_);
    
    // ========== Part 2: Coarse-grid setup ==========
    // TODO: Implement proper coarse-grid correction with Galerkin projection
    // For now, this is a placeholder - effectively making this a one-level method
    
    // Create simple coarse FE space
    int coarseOrder = 1;
    coarseFec_ = new mfem::H1_FECollection(coarseOrder, fineSpace_->GetMesh()->Dimension());
    coarseSpace_ = new mfem::ParFiniteElementSpace(fineSpace_->GetParMesh(), coarseFec_);
    
    coarseSize_ = coarseSpace_->GetTrueVSize();
    
    // Allocate temporary vectors
    coarseX_.SetSize(coarseSize_);
    coarseY_.SetSize(coarseSize_);
    tempFine_.SetSize(localSize_);
    
    // Build simple coarse operator by assembling on coarse space
    mfem::ConstantCoefficient one(1.0);
    mfem::ParBilinearForm* coarseForm = new mfem::ParBilinearForm(coarseSpace_);
    coarseForm->AddDomainIntegrator(new mfem::DiffusionIntegrator(one));
    coarseForm->Assemble();
    coarseForm->Finalize();
    coarseMatrix_ = coarseForm->ParallelAssemble();
    delete coarseForm;
    
    // Create AMG solver for coarse problem
    coarseSolver_ = new mfem::HypreBoomerAMG(*coarseMatrix_);
    dynamic_cast<mfem::HypreBoomerAMG*>(coarseSolver_)->SetPrintLevel(0);
    
    // Note: prolongation_ remains nullptr - proper implementation needed
}

TwoLevelSchwarz::~TwoLevelSchwarz()
{
    delete localSolver_;
    delete localMatrix_;
    
    delete coarseSolver_;
    delete coarseMatrix_;
    // Don't delete prolongation_ - it's owned by the FE space
    delete coarseSpace_;
    delete coarseFec_;
}

void TwoLevelSchwarz::Mult(const mfem::Vector& x, mfem::Vector& y) const
{
    // Two-level Schwarz: y = (Ψ A_H^{-1} Ψ^T + ∑ R_i^T A_i^{-1} R_i) x
    // TODO: Currently only local solves are implemented (one-level behavior)
    
    // ========== Part 1: Local subdomain solve ==========
    // y = ∑ R_i^T A_i^{-1} R_i x
    
    // Extract local portion of input vector
    const double* xData = x.GetData();
    double* localXData = localX_.GetData();
    for (int i = 0; i < localSize_; ++i)
    {
        localXData[i] = xData[i];
    }
    
    // Solve local system
    localY_ = 0.0;
    localSolver_->Mult(localX_, localY_);
    
    // Store local solution to output
    double* yData = y.GetData();
    const double* localYData = localY_.GetData();
    for (int i = 0; i < localSize_; ++i)
    {
        yData[i] = localYData[i];
    }
    
    // ========== Part 2: Coarse-grid correction (TODO) ==========
    // Currently disabled - proper restriction/prolongation needed
    // Future: Add y += Ψ A_H^{-1} Ψ^T x
}

void TwoLevelSchwarz::SetOperator(const mfem::Operator& op)
{
    // This preconditioner is constructed with a fixed operator
    // SetOperator is a no-op, provided for interface compliance
}

} // namespace hpcfem
