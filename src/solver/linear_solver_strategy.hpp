#ifndef MPFEM_LINEAR_SOLVER_STRATEGY_HPP
#define MPFEM_LINEAR_SOLVER_STRATEGY_HPP

#include "mpfem.hpp"

#include <memory>

namespace mpfem {

/**
 * @brief Abstract interface for linear solver strategies.
 * 
 * Error handling: errors are logged and the program exits immediately.
 */
class LinearSolverStrategy {
public:
    virtual ~LinearSolverStrategy() = default;

    /**
     * @brief Solve the linear system A x = b.
     */
    virtual void solve(mfem::SparseMatrix& matrix,
                       mfem::Vector& solution,
                       mfem::Vector& rhs) = 0;

    virtual void setMaxIterations(int maxIterations) = 0;
    virtual void setRelativeTolerance(double tolerance) = 0;
    virtual void setAbsoluteTolerance(double tolerance) = 0;
    virtual void setPrintLevel(int level) = 0;

    virtual int getNumIterations() const = 0;
    virtual double getFinalResidual() const = 0;
};

/**
 * @brief Factory function to create a solver strategy from configuration.
 */
std::unique_ptr<LinearSolverStrategy> createLinearSolver(
    const std::string& type,
    int maxIterations = 1000,
    double tolerance = 1e-10,
    int printLevel = 0);

} // namespace mpfem

#endif