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

    // Use direct solver
    // Note: MFEM's direct solver availability depends on build configuration
#ifdef MFEM_USE_SUITESPARSE
    mfem::UMFPackSolver umfSolver;
    umfSolver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;
    umfSolver.SetOperator(matrix);
    umfSolver.Mult(rhs, solution);
    numIterations_ = 1;
    finalResidual_ = 0.0;
#else
    Logger::log(LogLevel::Error, "Direct solver requested but UMFPACK is not available. Aborting.");
    std::exit(1);
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
