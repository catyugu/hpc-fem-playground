#ifndef MPFEM_CG_GS_SOLVER_HPP
#define MPFEM_CG_GS_SOLVER_HPP

#include "linear_solver_strategy.hpp"

namespace mpfem {

/**
 * @brief CG solver with Gauss-Seidel preconditioner.
 * 
 * This is the original solver implementation wrapped in the strategy pattern.
 * Suitable for small to medium problems where the Gauss-Seidel preconditioner
 * provides reasonable convergence.
 */
class CGGSSolver : public LinearSolverStrategy {
public:
    CGGSSolver();

    bool solve(mfem::SparseMatrix& matrix,
               mfem::Vector& solution,
               mfem::Vector& rhs,
               std::string& errorMessage) override;

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
