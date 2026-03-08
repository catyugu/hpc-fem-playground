#include "linear_solver_strategy.hpp"
#include "cg_gs_solver.hpp"

namespace mpfem {

std::unique_ptr<LinearSolverStrategy> createLinearSolver(
    const std::string& type,
    int maxIterations,
    double tolerance,
    int printLevel)
{
    std::unique_ptr<LinearSolverStrategy> solver;

    // All types now use the unified solver with appropriate settings
    auto cggs = std::make_unique<CGGSSolver>();
    cggs->setMaxIterations(maxIterations);
    cggs->setRelativeTolerance(tolerance);
    cggs->setPrintLevel(printLevel);
    solver = std::move(cggs);

    return solver;
}

} // namespace mpfem