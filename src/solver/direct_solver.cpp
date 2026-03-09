#include "direct_solver.hpp"
#include "mfem.hpp"
#include "logger.hpp"

#include <cmath>

namespace mpfem {

DirectSolver::DirectSolver()
    : maxIterations_(1000)
    , relativeTolerance_(1e-12)
    , absoluteTolerance_(0.0)
    , printLevel_(0)
    , numIterations_(0)
    , finalResidual_(0.0)
{
}

void DirectSolver::solve(mfem::SparseMatrix& matrix,
                         mfem::Vector& solution,
                         mfem::Vector& rhs)
{
    numIterations_ = 0;
    finalResidual_ = 0.0;

    // Validate RHS vector
    for (int i = 0; i < rhs.Size(); ++i) {
        Check(std::isfinite(rhs(i)),
              "DirectSolver: RHS vector contains non-finite value at index " +
              std::to_string(i) + " (value: " + std::to_string(rhs(i)) + ")");
    }

    // Try to use direct solver first
    // Note: MFEM's direct solver availability depends on build configuration
#ifdef MFEM_USE_SUITESPARSE
    // Use UMFPack if available
    mfem::UMFPackSolver umfSolver;
    umfSolver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;
    umfSolver.SetOperator(matrix);
    umfSolver.Mult(rhs, solution);
    numIterations_ = 1;
    finalResidual_ = 0.0;
#else
    // Fallback to iterative solver with diagonal scaling
    // Scale the system for better conditioning
    mfem::Vector diag(matrix.Size());
    matrix.GetDiag(diag);
    
    // Check for zeros on diagonal
    for (int i = 0; i < diag.Size(); ++i) {
        Check(std::abs(diag(i)) > 1e-15,
              "DirectSolver: Zero on diagonal at index " + std::to_string(i));
    }
    
    // Use MINRES with diagonal preconditioner
    mfem::Vector scaledRhs(rhs.Size());
    mfem::Vector scaledSol(solution.Size());
    
    for (int i = 0; i < rhs.Size(); ++i) {
        scaledRhs(i) = rhs(i) / std::sqrt(std::abs(diag(i)));
    }
    
    // Use simple iterative refinement
    solution = 0.0;
    mfem::Vector residual(rhs.Size());
    mfem::Vector correction(rhs.Size());
    
    mfem::GSSmoother gs(matrix);
    
    for (int iter = 0; iter < maxIterations_; ++iter) {
        // Compute residual: r = b - A*x
        matrix.Mult(solution, residual);
        for (int i = 0; i < rhs.Size(); ++i) {
            residual(i) = rhs(i) - residual(i);
        }
        
        double resNorm = residual.Norml2();
        if (iter == 0) {
            finalResidual_ = resNorm;
        }
        
        if (resNorm < relativeTolerance_ * finalResidual_ || resNorm < absoluteTolerance_) {
            numIterations_ = iter + 1;
            finalResidual_ = resNorm;
            break;
        }
        
        // Solve for correction using Gauss-Seidel
        correction = 0.0;
        gs.Mult(residual, correction);
        
        // Update solution
        for (int i = 0; i < solution.Size(); ++i) {
            solution(i) += correction(i);
        }
        
        numIterations_ = iter + 1;
        finalResidual_ = resNorm;
    }
#endif

    if (printLevel_ > 0) {
        Logger::log(LogLevel::Info, "DirectSolver: iterations=" + std::to_string(numIterations_) +
                ", residual=" + std::to_string(finalResidual_));
    }
}

void DirectSolver::setMaxIterations(int maxIterations)
{
    maxIterations_ = maxIterations;
}

void DirectSolver::setRelativeTolerance(double tolerance)
{
    relativeTolerance_ = tolerance;
}

void DirectSolver::setAbsoluteTolerance(double tolerance)
{
    absoluteTolerance_ = tolerance;
}

void DirectSolver::setPrintLevel(int level)
{
    printLevel_ = level;
}

int DirectSolver::getNumIterations() const
{
    return numIterations_;
}

double DirectSolver::getFinalResidual() const
{
    return finalResidual_;
}

} // namespace mpfem
