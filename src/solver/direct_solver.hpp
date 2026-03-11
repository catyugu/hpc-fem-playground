#ifndef MPFEM_DIRECT_SOLVER_HPP
#define MPFEM_DIRECT_SOLVER_HPP

#include "linear_solver_strategy.hpp"

namespace mpfem {

/**
 * @brief Direct solver using UMFPACK (SuiteSparse).
 */
class DirectSolver : public LinearSolverStrategy {
public:
    void solve(mfem::SparseMatrix& matrix,
               mfem::Vector& solution,
               mfem::Vector& rhs) override;

    // These are ignored for direct solvers
    void setMaxIterations(int) override {}
    void setRelativeTolerance(double) override {}
    void setAbsoluteTolerance(double) override {}
    void setPrintLevel(int level) override { printLevel_ = level; }

    int getNumIterations() const override { return 1; }
    double getFinalResidual() const override { return 0.0; }

private:
    int printLevel_ = 0;
};

} // namespace mpfem

#endif