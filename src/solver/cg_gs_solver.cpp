#include "cg_gs_solver.hpp"
#include "mpfem_types.hpp"
#include "logger.hpp"

namespace mpfem {

CGGSSolver::CGGSSolver()
    : maxIterations_(1000)
    , relativeTolerance_(1e-10)
    , absoluteTolerance_(0.0)
    , printLevel_(0)
    , numIterations_(0)
    , finalResidual_(0.0)
{
}

void CGGSSolver::solve(FemMatrix& matrix,
                       FemVector& solution,
                       FemVector& rhs)
{
    numIterations_ = 0;
    finalResidual_ = 0.0;

#ifdef MFEM_USE_MPI
    // Parallel version with HYPRE AMG preconditioner
    mfem::HypreParMatrix* hypreMat = dynamic_cast<mfem::HypreParMatrix*>(&matrix);
    Check(hypreMat != nullptr, "Matrix must be HypreParMatrix for parallel solver");

    mfem::HypreBoomerAMG amg(*hypreMat);
    amg.SetPrintLevel(printLevel_);

    // Use MFEM's CGSolver which works with any Operator
    mfem::CGSolver cg(hypreMat->GetComm());
    cg.SetOperator(*hypreMat);
    cg.SetPreconditioner(amg);
    cg.SetMaxIter(maxIterations_);
    cg.SetRelTol(relativeTolerance_);
    cg.SetAbsTol(absoluteTolerance_);
    cg.SetPrintLevel(printLevel_);
    cg.Mult(rhs, solution);

    numIterations_ = cg.GetNumIterations();
    finalResidual_ = cg.GetFinalNorm();

    Check(cg.GetConverged(), 
          "CG solver did not converge after " + std::to_string(numIterations_) 
          + " iterations. Final residual: " + std::to_string(finalResidual_));
#else
    // Serial version with CG + Gauss-Seidel
    mfem::GSSmoother preconditioner(matrix);
    mfem::CGSolver cgSolver;
    cgSolver.SetPreconditioner(preconditioner);
    cgSolver.SetOperator(matrix);
    cgSolver.SetMaxIter(maxIterations_);
    cgSolver.SetRelTol(relativeTolerance_);
    cgSolver.SetAbsTol(absoluteTolerance_);
    cgSolver.SetPrintLevel(printLevel_);
    cgSolver.Mult(rhs, solution);

    numIterations_ = cgSolver.GetNumIterations();
    finalResidual_ = cgSolver.GetFinalNorm();

    Check(cgSolver.GetConverged(), 
          "CG solver did not converge after " + std::to_string(numIterations_) 
          + " iterations. Final residual: " + std::to_string(finalResidual_));
#endif
}

void CGGSSolver::setMaxIterations(int maxIterations)
{
    maxIterations_ = maxIterations;
}

void CGGSSolver::setRelativeTolerance(double tolerance)
{
    relativeTolerance_ = tolerance;
}

void CGGSSolver::setAbsoluteTolerance(double tolerance)
{
    absoluteTolerance_ = tolerance;
}

void CGGSSolver::setPrintLevel(int level)
{
    printLevel_ = level;
}

int CGGSSolver::getNumIterations() const
{
    return numIterations_;
}

double CGGSSolver::getFinalResidual() const
{
    return finalResidual_;
}

} // namespace mpfem
