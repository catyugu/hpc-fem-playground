#include "linear_solver_strategy.hpp"
#include "cg_gs_solver.hpp"
#include "direct_solver.hpp"
#include "logger.hpp"

namespace mpfem {

std::unique_ptr<LinearSolverStrategy> createLinearSolver(
    const std::string& type,
    int maxIterations,
    double tolerance,
    int printLevel)
{
    std::unique_ptr<LinearSolverStrategy> solver;

    if (type == "direct" || type == "umfpack") {
        // Direct solver
#ifdef MFEM_USE_SUITESPARSE
        auto directSolver = std::make_unique<DirectSolver>();
        directSolver->setMaxIterations(maxIterations);
        directSolver->setRelativeTolerance(tolerance);
        directSolver->setPrintLevel(printLevel);
        solver = std::move(directSolver);
        Logger::log(LogLevel::Info, "Using direct solver");
#else
        Logger::log(LogLevel::Error, "Direct solver requested but UMFPACK is not available. Aborting.");
        std::exit(1);
#endif
    } else if (type == "cg_gs") {
        // CG with Gauss-Seidel preconditioner
        auto cggs = std::make_unique<CGGSSolver>();
        cggs->setMaxIterations(maxIterations);
        cggs->setRelativeTolerance(tolerance);
        cggs->setPrintLevel(printLevel);
        solver = std::move(cggs);
        Logger::log(LogLevel::Info, "Using CG solver with Gauss-Seidel preconditioner");
    } else{
        Logger::log(LogLevel::Error, "Unknown linear solver type: " + type + ". Aborting.");
        std::exit(1);
    }

    return solver;
}

} // namespace mpfem