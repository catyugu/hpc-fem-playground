/**
 * @file solvers/solver_hypre_amg.cpp
 * @brief Implementation of HYPRE BoomerAMG solver
 */

#include "hpcfem/solvers/solver_hypre_amg.hpp"

namespace hpcfem
{

HypreAmgSolver::HypreAmgSolver(double relTol, int maxIter, int printLevel)
    : relTol_(relTol), maxIter_(maxIter), printLevel_(printLevel), numIterations_(0), finalNorm_(0.0)
{
}

#ifdef MFEM_USE_MPI
void HypreAmgSolver::solve(const mfem::HypreParMatrix& A,
                           const mfem::Vector& b,
                           mfem::Vector& x)
{
    mfem::HypreParMatrix* A_copy = const_cast<mfem::HypreParMatrix*>(&A);
    mfem::HypreBoomerAMG prec(*A_copy);
    prec.SetPrintLevel(0);

    mfem::CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(relTol_);
    cg.SetMaxIter(maxIter_);
    cg.SetPrintLevel(printLevel_);
    cg.SetPreconditioner(prec);
    cg.SetOperator(*A_copy);

    cg.Mult(b, x);

    numIterations_ = cg.GetNumIterations();
    finalNorm_ = cg.GetFinalNorm();
}
#else
void HypreAmgSolver::solve(const mfem::SparseMatrix& A,
                           const mfem::Vector& b,
                           mfem::Vector& x)
{
    mfem::SparseMatrix* A_copy = const_cast<mfem::SparseMatrix*>(&A);
    mfem::GSSmoother prec(*A_copy);

    mfem::CGSolver cg;
    cg.SetRelTol(relTol_);
    cg.SetMaxIter(maxIter_);
    cg.SetPrintLevel(printLevel_);
    cg.SetPreconditioner(prec);
    cg.SetOperator(*A_copy);

    cg.Mult(b, x);

    numIterations_ = cg.GetNumIterations();
    finalNorm_ = cg.GetFinalNorm();
}
#endif

} // namespace hpcfem
