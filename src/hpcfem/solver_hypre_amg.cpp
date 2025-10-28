/**
 * @file solver_hypre_amg.cpp
 * @brief Implementation of HYPRE BoomerAMG solver
 * 
 * Phase: 2, Step: 2.2
 */

#include "solver_hypre_amg.hpp"

namespace hpcfem
{

HypreAmgSolver::HypreAmgSolver(double relTol, int maxIter, int printLevel)
    : relTol_(relTol)
    , maxIter_(maxIter)
    , printLevel_(printLevel)
    , numIterations_(0)
    , finalNorm_(0.0)
{
}

#ifdef MFEM_USE_MPI
void HypreAmgSolver::solve(const mfem::HypreParMatrix& A,
                           const mfem::Vector& b,
                           mfem::Vector& x)
{
    // Create a mutable copy of the matrix for solver operations
    mfem::HypreParMatrix* A_copy = const_cast<mfem::HypreParMatrix*>(&A);
    
    // Create AMG preconditioner
    mfem::HypreBoomerAMG prec(*A_copy);
    prec.SetPrintLevel(0);  // Suppress AMG setup output
    
    // Create CG solver
    mfem::CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(relTol_);
    cg.SetMaxIter(maxIter_);
    cg.SetPrintLevel(printLevel_);
    cg.SetPreconditioner(prec);
    cg.SetOperator(*A_copy);
    
    // Solve the system
    cg.Mult(b, x);
    
    // Store solver statistics
    numIterations_ = cg.GetNumIterations();
    finalNorm_ = cg.GetFinalNorm();
}
#else
void HypreAmgSolver::solve(const mfem::SparseMatrix& A,
                           const mfem::Vector& b,
                           mfem::Vector& x)
{
    // Create a mutable copy of the matrix for solver operations
    mfem::SparseMatrix* A_copy = const_cast<mfem::SparseMatrix*>(&A);
    
    // Create Gauss-Seidel preconditioner (no HYPRE in serial mode)
    mfem::GSSmoother prec(*A_copy);
    
    // Create CG solver
    mfem::CGSolver cg;
    cg.SetRelTol(relTol_);
    cg.SetMaxIter(maxIter_);
    cg.SetPrintLevel(printLevel_);
    cg.SetPreconditioner(prec);
    cg.SetOperator(*A_copy);
    
    // Solve the system
    cg.Mult(b, x);
    
    // Store solver statistics
    numIterations_ = cg.GetNumIterations();
    finalNorm_ = cg.GetFinalNorm();
}
#endif

} // namespace hpcfem
