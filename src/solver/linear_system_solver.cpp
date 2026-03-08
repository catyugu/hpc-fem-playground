#include "linear_system_solver.hpp"

namespace mpfem {

bool LinearSystemSolver::solve(mfem::SparseMatrix &matrix,
                               mfem::Vector &rhs,
                               mfem::Vector &solution,
                               std::string &errorMessage)
{
    errorMessage.clear();

    mfem::GSSmoother preconditioner(matrix);
    mfem::PCG(matrix,
              preconditioner,
              rhs,
              solution,
              0,
              500,
              1e-12,
              0.0);

    if (solution.Size() != rhs.Size()) {
        errorMessage = "LinearSystemSolver produced invalid solution size";
        return false;
    }

    return true;
}

} // namespace mpfem
