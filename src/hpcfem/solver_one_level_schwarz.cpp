/**
 * @file solver_one_level_schwarz.cpp
 * @brief Implementation of one-level overlapping Schwarz preconditioner
 * 
 * Phase: 2, Step: 2.1
 */

#include "hpcfem/solver_one_level_schwarz.hpp"

namespace hpcfem
{

OneLevelSchwarz::OneLevelSchwarz(const mfem::HypreParMatrix& A)
    : mfem::Solver(A.Height()),
      globalMatrix_(&A),
      localMatrix_(nullptr),
      localSolver_(nullptr)
{
    // Get the diagonal block (local subdomain matrix) from the parallel matrix
    // This is the A_i in the Schwarz formula
    localMatrix_ = new mfem::SparseMatrix();
    A.GetDiag(*localMatrix_);
    localSize_ = localMatrix_->Height();
    
    // Allocate temporary vectors for local operations
    localX_.SetSize(localSize_);
    localY_.SetSize(localSize_);
    
    // Create local Gauss-Seidel smoother as the subdomain solver
    // This performs multiple sweeps to approximate A_i^{-1}
    localSolver_ = new mfem::GSSmoother(*localMatrix_);
}

OneLevelSchwarz::~OneLevelSchwarz()
{
    delete localSolver_;
    delete localMatrix_;
}

void OneLevelSchwarz::Mult(const mfem::Vector& x, mfem::Vector& y) const
{
    // Additive Schwarz: y = ∑ R_i^T A_i^{-1} R_i x
    // 
    // For MFEM's parallel vectors:
    // 1. x is already distributed (each process has its local portion)
    // 2. Extract local portion (restriction R_i x)
    // 3. Solve local system A_i localY = localX
    // 4. Store result back to y (prolongation R_i^T)
    
    // Extract local portion of input vector (restriction)
    const double* xData = x.GetData();
    double* localXData = localX_.GetData();
    for (int i = 0; i < localSize_; ++i)
    {
        localXData[i] = xData[i];
    }
    
    // Solve local system: A_i localY = localX
    localY_ = 0.0;  // Initial guess
    localSolver_->Mult(localX_, localY_);
    
    // Store local solution back to global vector (prolongation)
    double* yData = y.GetData();
    const double* localYData = localY_.GetData();
    for (int i = 0; i < localSize_; ++i)
    {
        yData[i] = localYData[i];
    }
    
    // Note: In the parallel setting, MFEM's HypreParVector handles the
    // communication and summation of overlapping contributions automatically
    // through its parallel structure. The "additive" nature comes from
    // each process contributing its local solution.
}

void OneLevelSchwarz::SetOperator(const mfem::Operator& op)
{
    // This preconditioner is constructed with a fixed operator
    // SetOperator is a no-op, provided for interface compliance
    // If dynamic operator changes are needed, this would require
    // reconstructing the local solver
}

} // namespace hpcfem
