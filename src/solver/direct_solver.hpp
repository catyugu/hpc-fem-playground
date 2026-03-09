#ifndef MPFEM_DIRECT_SOLVER_HPP
#define MPFEM_DIRECT_SOLVER_HPP

#include "linear_solver_strategy.hpp"

namespace mpfem {

/**
 * @brief Direct solver using MFEM's built-in capabilities.
 * Falls back to iterative solver if direct solver is not available.
 */
class DirectSolver : public LinearSolverStrategy {
public:
    DirectSolver();

    void solve(mfem::SparseMatrix& matrix,
               mfem::Vector& solution,
               mfem::Vector& rhs) override;

    void setMaxIterations(int maxIterations) override;
    void setRelativeTolerance(double tolerance) override;
    void setAbsoluteTolerance(double tolerance) override;
    void setPrintLevel(int level) override;

    int getNumIterations() const override;
    double getFinalResidual() const override;

private:
    int maxIterations_;
    double relativeTolerance_;
    double absoluteTolerance_;
    int printLevel_;
    mutable int numIterations_;
    mutable double finalResidual_;
};

} // namespace mpfem

#endif
