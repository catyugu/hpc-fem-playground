/**
 * @file solver_ddm_schwarz.cpp
 * @brief Implementation of DDM Schwarz solver
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
      omega_(omega),
      numIterations_(0),
      finalNorm_(0.0)
{
    // Note: Optimal omega for Richardson iteration is typically 0.5-0.7
    // omega = 1.0 can be too aggressive and cause slow convergence
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
    
    // Setup local subdomain solver (AMG preconditioner)
    // CRITICAL: Build the AMG hierarchy ONCE, not every iteration!
    mfem::HypreBoomerAMG localSolver;
    localSolver.SetPrintLevel(0);
    localSolver.SetOperator(A);  // Build AMG hierarchy once here
    
    // Create temporary vectors
    mfem::Vector r(b.Size());  // residual
    mfem::Vector z(b.Size());  // preconditioned residual
    
    // Compute initial residual: r = b - A*x
    A.Mult(x, r);
    subtract(b, r, r);
    
    double initialNorm = mfem::InnerProduct(A.GetComm(), r, r);
    initialNorm = std::sqrt(initialNorm);
    
    if (printLevel_ > 0 && myid == 0)
    {
        std::cout << "DDM Schwarz: Initial residual norm = " << initialNorm << std::endl;
    }
    
    double normThreshold = relTol_ * initialNorm;
    double currentNorm = initialNorm;
    
    // Richardson iteration with Schwarz preconditioning
    for (int iter = 0; iter < maxIter_; ++iter)
    {
        // Apply local subdomain solver: z = M^{-1} * r
        // The AMG hierarchy is already built, just apply it
        z = 0.0;
        localSolver.Mult(r, z);
        
        // Update solution: x = x + omega * z
        x.Add(omega_, z);
        
        // Compute new residual: r = b - A*x
        A.Mult(x, r);
        subtract(b, r, r);
        
        // Compute global L2 norm of residual
        currentNorm = mfem::InnerProduct(A.GetComm(), r, r);
        currentNorm = std::sqrt(currentNorm);
        
        if (printLevel_ > 1 && myid == 0)
        {
            std::cout << "  Iteration " << iter + 1 
                      << ": residual norm = " << currentNorm << std::endl;
        }
        
        // Check convergence
        if (currentNorm < normThreshold)
        {
            numIterations_ = iter + 1;
            finalNorm_ = currentNorm;
            
            if (printLevel_ > 0 && myid == 0)
            {
                std::cout << "DDM Schwarz converged in " << numIterations_ 
                          << " iterations (final norm: " << finalNorm_ << ")" << std::endl;
            }
            return;
        }
    }
    
    // Did not converge
    numIterations_ = maxIter_;
    finalNorm_ = currentNorm;
    
    if (printLevel_ > 0 && myid == 0)
    {
        std::cout << "DDM Schwarz: Maximum iterations reached. Final norm: " 
                  << finalNorm_ << std::endl;
    }
}
#else
void DdmSchwarzSolver::solve(const mfem::SparseMatrix& A,
                             const mfem::Vector& b,
                             mfem::Vector& x)
{
    // Serial version: Use GS smoother as local solver
    // Build the smoother once, not every iteration
    mfem::GSSmoother localSolver(A);
    
    mfem::Vector r(b.Size());
    mfem::Vector z(b.Size());
    
    // Compute initial residual
    A.Mult(x, r);
    subtract(b, r, r);
    
    double initialNorm = r.Norml2();
    
    if (printLevel_ > 0)
    {
        std::cout << "DDM Schwarz (serial): Initial residual norm = " << initialNorm << std::endl;
    }
    
    double normThreshold = relTol_ * initialNorm;
    double currentNorm = initialNorm;
    
    // Richardson iteration
    for (int iter = 0; iter < maxIter_; ++iter)
    {
        // Apply smoother (already initialized)
        z = 0.0;
        localSolver.Mult(r, z);
        
        x.Add(omega_, z);
        
        A.Mult(x, r);
        subtract(b, r, r);
        
        currentNorm = r.Norml2();
        
        if (printLevel_ > 1)
        {
            std::cout << "  Iteration " << iter + 1 
                      << ": residual norm = " << currentNorm << std::endl;
        }
        
        if (currentNorm < normThreshold)
        {
            numIterations_ = iter + 1;
            finalNorm_ = currentNorm;
            
            if (printLevel_ > 0)
            {
                std::cout << "DDM Schwarz (serial) converged in " << numIterations_ 
                          << " iterations" << std::endl;
            }
            return;
        }
    }
    
    numIterations_ = maxIter_;
    finalNorm_ = currentNorm;
    
    if (printLevel_ > 0)
    {
        std::cout << "DDM Schwarz (serial): Maximum iterations reached" << std::endl;
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
