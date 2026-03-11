#include "pardiso_solver.hpp"
#include "logger.hpp"

#ifdef MFEM_USE_MKL_PARDISO
#include "mfem.hpp"
#endif

namespace mpfem {

void PardisoSolver::solve(mfem::SparseMatrix& matrix,
                          mfem::Vector& solution,
                          mfem::Vector& rhs)
{
#ifdef MFEM_USE_MKL_PARDISO
    mfem::PardisoSolver pardiso;
    pardiso.SetPrintLevel(printLevel_);
    pardiso.SetMatrixType(mfem::PardisoSolver::MatType::REAL_NONSYMMETRIC);
    pardiso.SetOperator(matrix);
    pardiso.Mult(rhs, solution);
#else
    Logger::log(LogLevel::Error, "PARDISO solver requested but MKL PARDISO is not available. Aborting.");
    std::exit(1);
#endif
}

} // namespace mpfem
