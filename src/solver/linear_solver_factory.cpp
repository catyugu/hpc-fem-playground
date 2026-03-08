#include "linear_solver_strategy.hpp"
#include "cg_gs_solver.hpp"
#include "hypre_boomeramg_solver.hpp"

#include <stdexcept>

namespace mpfem {

std::unique_ptr<LinearSolverStrategy> createLinearSolver(
    const std::string& type,
    int maxIterations,
    double tolerance,
    int printLevel)
{
    std::unique_ptr<LinearSolverStrategy> solver;

    if (type == "cg_gs") {
        auto cggs = std::make_unique<CGGSSolver>();
        cggs->setMaxIterations(maxIterations);
        cggs->setRelativeTolerance(tolerance);
        cggs->setPrintLevel(printLevel);
        solver = std::move(cggs);
    }
    else if (type == "hypre_boomeramg" || type == "hypre_amg") {
        bool useOuterIteration = (type == "hypre_boomeramg");
        auto hypreSolver = std::make_unique<HYPREBoomerAMGSolver>(useOuterIteration);
        hypreSolver->setMaxIterations(maxIterations);
        hypreSolver->setRelativeTolerance(tolerance);
        hypreSolver->setPrintLevel(printLevel);
        solver = std::move(hypreSolver);
    }
    else if (type == "gmres") {
        // Use HYPRE AMG preconditioner with GMRES (good for non-symmetric problems)
        auto hypreSolver = std::make_unique<HYPREBoomerAMGSolver>(true);
        hypreSolver->setMaxIterations(maxIterations);
        hypreSolver->setRelativeTolerance(tolerance);
        hypreSolver->setPrintLevel(printLevel);
        solver = std::move(hypreSolver);
    }
    else {
        // Default to CG+GS for unknown types
        auto cggs = std::make_unique<CGGSSolver>();
        cggs->setMaxIterations(maxIterations);
        cggs->setRelativeTolerance(tolerance);
        cggs->setPrintLevel(printLevel);
        solver = std::move(cggs);
    }

    return solver;
}

} // namespace mpfem
