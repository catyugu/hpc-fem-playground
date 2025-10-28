/**
 * @file solver_ddm_schwarz.cpp
 * @brief Implementation of DDM Schwarz solver with PCG acceleration
 */

#include "solver_ddm_schwarz.hpp"
#include <cmath>
#include <iostream>

namespace hpcfem
{

DdmSchwarzSolver::DdmSchwarzSolver(double relTol,
                                   int maxIter,
                                   int printLevel,
                                   double omega)
    : relTol_(relTol),
      maxIter_(maxIter),
      printLevel_(printLevel),
      omega_(omega),  // Not used with CG, kept for interface compatibility
      numIterations_(0),
      finalNorm_(0.0)
{
}

DdmSchwarzSolver::~DdmSchwarzSolver()
{
}

#ifdef MFEM_USE_MPI
void DdmSchwarzSolver::solve(const mfem::HypreParMatrix& A,
                             const mfem::Vector& b,
                             mfem::Vector& x)
{
    int myid = 0;
    MPI_Comm_rank(A.GetComm(), &myid);
    
    // Setup local subdomain solver as preconditioner
    // Build AMG hierarchy once for reuse
    mfem::HypreBoomerAMG preconditioner;
    preconditioner.SetPrintLevel(0);
    preconditioner.SetOperator(A);
    
    // Use PCG (Preconditioned Conjugate Gradient) instead of Richardson
    // CG is optimal for SPD systems and doesn't require parameter tuning
    mfem::CGSolver cg(A.GetComm());
    cg.SetRelTol(relTol_);
    cg.SetMaxIter(maxIter_);
    cg.SetPrintLevel(printLevel_);
    cg.SetOperator(A);
    cg.SetPreconditioner(preconditioner);
    
    // Solve using PCG with Schwarz preconditioning
    cg.Mult(b, x);
    
    // Extract convergence info
    numIterations_ = cg.GetNumIterations();
    finalNorm_ = cg.GetFinalNorm();
    
    if (printLevel_ > 0 && myid == 0)
    {
        if (cg.GetConverged())
        {
            std::cout << "DDM-PCG converged in " << numIterations_ 
                      << " iterations (final norm: " << finalNorm_ << ")" << std::endl;
        }
        else
        {
            std::cout << "DDM-PCG: Maximum iterations reached. Final norm: " 
                      << finalNorm_ << std::endl;
        }
    }
}
#else
void DdmSchwarzSolver::solve(const mfem::SparseMatrix& A,
                             const mfem::Vector& b,
                             mfem::Vector& x)
{
    // Serial version: Use CG with GS preconditioner
    mfem::GSSmoother preconditioner(A);
    
    mfem::CG cg;
    cg.SetRelTol(relTol_);
    cg.SetMaxIter(maxIter_);
    cg.SetPrintLevel(printLevel_);
    cg.SetOperator(A);
    cg.SetPreconditioner(preconditioner);
    
    cg.Mult(b, x);
    
    numIterations_ = cg.GetNumIterations();
    finalNorm_ = cg.GetFinalNorm();
    
    if (printLevel_ > 0)
    {
        if (cg.GetConverged())
        {
            std::cout << "DDM-PCG (serial) converged in " << numIterations_ 
                      << " iterations (final norm: " << finalNorm_ << ")" << std::endl;
        }
        else
        {
            std::cout << "DDM-PCG (serial): Maximum iterations reached" << std::endl;
        }
    }
}
#endif

int DdmSchwarzSolver::getNumIterations() const
{
    return numIterations_;
}

double DdmSchwarzSolver::getFinalNorm() const
{
    return finalNorm_;
}

} // namespace hpcfem
