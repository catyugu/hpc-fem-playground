#include "cg_gs_solver.hpp"
#include "mfem.hpp"
#include "logger.hpp"

#include <cmath>

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

void CGGSSolver::solve(mfem::SparseMatrix& matrix,
                       mfem::Vector& solution,
                       mfem::Vector& rhs)
{
    numIterations_ = 0;
    finalResidual_ = 0.0;

    // Validate RHS vector
    for (int i = 0; i < rhs.Size(); ++i) {
        Check(std::isfinite(rhs(i)),
              "CGGSSolver: RHS vector contains non-finite value at index " +
              std::to_string(i) + " (value: " + std::to_string(rhs(i)) + ")");
    }

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
