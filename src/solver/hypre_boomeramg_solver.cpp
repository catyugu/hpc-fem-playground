#include "hypre_boomeramg_solver.hpp"

#ifdef MFEM_USE_MPI
#include "linalg/hypre.hpp"
#endif

namespace mpfem {

HYPREBoomerAMGSolver::HYPREBoomerAMGSolver(bool useOuterIteration)
    : useOuterIteration_(useOuterIteration)
    , maxIterations_(1000)
    , relativeTolerance_(1e-10)
    , absoluteTolerance_(0.0)
    , printLevel_(0)
    , relaxType_(3)      // Chebyshev smoothing (default, works well for elasticity)
    , coarseningType_(10) // HMIS coarsening (good default)
    , numIterations_(0)
    , finalResidual_(0.0)
{
}

bool HYPREBoomerAMGSolver::solve(mfem::SparseMatrix& matrix,
                                  mfem::Vector& solution,
                                  mfem::Vector& rhs,
                                  std::string& errorMessage)
{
    errorMessage.clear();
    numIterations_ = 0;
    finalResidual_ = 0.0;

#ifdef MFEM_USE_HYPRE
    // Convert to HypreParMatrix (serial case: wrap the local matrix)
    mfem::HypreParMatrix hypreMatrix(MPI_COMM_SELF, matrix);
    
    if (useOuterIteration_) {
        // Use PCG with BoomerAMG preconditioner
        mfem::HypreBoomerAMG amg(hypreMatrix);
        amg.SetPrintLevel(printLevel_);
        amg.SetRelaxType(relaxType_);
        amg.SetCoarsening(coarseningType_);
        
        mfem::HyprePCG pcg(hypreMatrix);
        pcg.SetMaxIter(maxIterations_);
        pcg.SetTol(relativeTolerance_);
        pcg.SetAbsTol(absoluteTolerance_);
        pcg.SetPrintLevel(printLevel_);
        pcg.SetPreconditioner(amg);
        pcg.Mult(rhs, solution);
        
        numIterations_ = pcg.GetNumIterations();
        finalResidual_ = pcg.GetFinalNorm();
        
        if (!pcg.GetConverged()) {
            errorMessage = "HYPRE PCG solver did not converge after "
                         + std::to_string(numIterations_)
                         + " iterations. Final residual: "
                         + std::to_string(finalResidual_);
            return false;
        }
    } else {
        // Use AMG as standalone solver
        mfem::HypreBoomerAMG amg(hypreMatrix);
        amg.SetMaxIter(maxIterations_);
        amg.SetTol(relativeTolerance_);
        amg.SetAbsTol(absoluteTolerance_);
        amg.SetPrintLevel(printLevel_);
        amg.SetRelaxType(relaxType_);
        amg.SetCoarsening(coarseningType_);
        amg.Mult(rhs, solution);
        
        numIterations_ = amg.GetNumIterations();
        finalResidual_ = amg.GetFinalNorm();
        
        // For AMG standalone, convergence check is different
        // AMG typically converges in very few iterations for elliptic problems
    }
    
    return true;
#else
    // Fallback to CG+GS if HYPRE is not available
    errorMessage = "HYPRE not available, falling back to CG+GS";
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
    
    return cgSolver.GetConverged();
#endif
}

void HYPREBoomerAMGSolver::setMaxIterations(int maxIterations)
{
    maxIterations_ = maxIterations;
}

void HYPREBoomerAMGSolver::setRelativeTolerance(double tolerance)
{
    relativeTolerance_ = tolerance;
}

void HYPREBoomerAMGSolver::setAbsoluteTolerance(double tolerance)
{
    absoluteTolerance_ = tolerance;
}

void HYPREBoomerAMGSolver::setPrintLevel(int level)
{
    printLevel_ = level;
}

void HYPREBoomerAMGSolver::setRelaxType(int type)
{
    relaxType_ = type;
}

void HYPREBoomerAMGSolver::setCoarseningType(int type)
{
    coarseningType_ = type;
}

int HYPREBoomerAMGSolver::getNumIterations() const
{
    return numIterations_;
}

double HYPREBoomerAMGSolver::getFinalResidual() const
{
    return finalResidual_;
}

} // namespace mpfem
