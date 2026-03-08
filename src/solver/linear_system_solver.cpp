#include "linear_system_solver.hpp"

namespace mpfem {

bool LinearSystemSolver::solve(mfem::SparseMatrix &matrix,
                               mfem::Vector &rhs,
                               mfem::Vector &solution,
                               int maxIterations,
                               double relativeTolerance,
                               std::string &errorMessage)
{
    errorMessage.clear();

    mfem::GSSmoother preconditioner(matrix);
    mfem::CGSolver cgSolver;
    cgSolver.SetOperator(matrix);
    cgSolver.SetPreconditioner(preconditioner);
    cgSolver.SetMaxIter(maxIterations > 0 ? maxIterations : 500);
    cgSolver.SetRelTol(relativeTolerance > 0.0 ? relativeTolerance : 1e-10);
    cgSolver.SetAbsTol(0.0);
    cgSolver.SetPrintLevel(0);
    cgSolver.Mult(rhs, solution);

    if (solution.Size() != rhs.Size()) {
        errorMessage = "LinearSystemSolver produced invalid solution size";
        return false;
    }

    if (!cgSolver.GetConverged()) {
        errorMessage = "LinearSystemSolver failed to converge in CG";
        return false;
    }

    return true;
}

} // namespace mpfem
