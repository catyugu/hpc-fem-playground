#ifndef MPFEM_PARDISO_SOLVER_HPP
#define MPFEM_PARDISO_SOLVER_HPP

#include "linear_solver_strategy.hpp"

namespace mpfem {

/**
 * @brief Direct solver using Intel MKL PARDISO.
 * 
 * PARDISO is typically faster than UMFPACK for larger problems and supports
 * multi-threaded parallelism within the solver.
 */
class PardisoSolver : public LinearSolverStrategy {
public:
    void solve(mfem::SparseMatrix& matrix,
               mfem::Vector& solution,
               mfem::Vector& rhs) override;

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
