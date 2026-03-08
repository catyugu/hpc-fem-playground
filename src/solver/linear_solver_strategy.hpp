#ifndef MPFEM_LINEAR_SOLVER_STRATEGY_HPP
#define MPFEM_LINEAR_SOLVER_STRATEGY_HPP

#include "mfem.hpp"

#include <memory>
#include <string>

namespace mpfem {

/**
 * @brief Abstract interface for linear solver strategies.
 * 
 * This class defines the interface for different linear solver implementations,
 * allowing the physics solvers to use different solving strategies (CG, GMRES,
 * HYPRE AMG, etc.) through a unified interface.
 */
class LinearSolverStrategy {
public:
    virtual ~LinearSolverStrategy() = default;

    /**
     * @brief Solve the linear system A x = b.
     * @param matrix System matrix (sparse).
     * @param solution Solution vector (input: initial guess, output: solution).
     * @param rhs Right-hand side vector.
     * @param errorMessage Error message on failure.
     * @return True if solved successfully.
     */
    virtual bool solve(mfem::SparseMatrix& matrix,
                       mfem::Vector& solution,
                       mfem::Vector& rhs,
                       std::string& errorMessage) = 0;

    /**
     * @brief Set the maximum number of iterations.
     */
    virtual void setMaxIterations(int maxIterations) = 0;

    /**
     * @brief Set the relative tolerance for convergence.
     */
    virtual void setRelativeTolerance(double tolerance) = 0;

    /**
     * @brief Set the absolute tolerance for convergence.
     */
    virtual void setAbsoluteTolerance(double tolerance) = 0;

    /**
     * @brief Set the print level for solver output.
     */
    virtual void setPrintLevel(int level) = 0;

    /**
     * @brief Get the number of iterations from the last solve.
     */
    virtual int getNumIterations() const = 0;

    /**
     * @brief Get the final residual from the last solve.
     */
    virtual double getFinalResidual() const = 0;
};

/**
 * @brief Factory function to create a solver strategy from configuration.
 * @param type Solver type: "cg_gs", "hypre_boomeramg", "hypre_amg", "gmres".
 * @param maxIterations Maximum iterations.
 * @param tolerance Relative tolerance.
 * @param printLevel Print level.
 * @return Unique pointer to solver strategy.
 */
std::unique_ptr<LinearSolverStrategy> createLinearSolver(
    const std::string& type,
    int maxIterations = 1000,
    double tolerance = 1e-10,
    int printLevel = 0);

} // namespace mpfem

#endif
