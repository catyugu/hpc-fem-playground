#include "direct_solver.hpp"
#include "mfem.hpp"
#include "logger.hpp"

#include <cmath>

namespace mpfem {

void DirectSolver::solve(mfem::SparseMatrix& matrix,
                         mfem::Vector& solution,
                         mfem::Vector& rhs)
{
    // Validate RHS vector
    for (int i = 0; i < rhs.Size(); ++i) {
        Check(std::isfinite(rhs(i)),
              "DirectSolver: RHS contains non-finite value at index " +
              std::to_string(i) + " (value: " + std::to_string(rhs(i)) + ")");
    }

#ifdef MFEM_USE_SUITESPARSE
    mfem::UMFPackSolver umfSolver;
    umfSolver.Control[UMFPACK_ORDERING] = UMFPACK_ORDERING_BEST;
    umfSolver.SetOperator(matrix);
    umfSolver.Mult(rhs, solution);
#else
    Logger::log(LogLevel::Error, "Direct solver requires UMFPACK. Aborting.");
    std::exit(1);
#endif

    if (printLevel_ > 0) {
        Logger::log(LogLevel::Info, "DirectSolver: solved");
    }
}

} // namespace mpfem