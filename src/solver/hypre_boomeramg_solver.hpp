#ifndef MPFEM_HYPRE_BOOMERAMG_SOLVER_HPP
#define MPFEM_HYPRE_BOOMERAMG_SOLVER_HPP

#include "linear_solver_strategy.hpp"

namespace mpfem {

/**
 * @brief HYPRE BoomerAMG solver with optional CG or GMRES outer iteration.
 * 
 * Uses algebraic multigrid (AMG) from HYPRE as a preconditioner for CG or GMRES,
 * or as a standalone solver. This is highly efficient for elliptic PDEs and
 * scales well for large problems.
 */
class HYPREBoomerAMGSolver : public LinearSolverStrategy {
public:
    /**
     * @brief Constructor.
     * @param useOuterIteration If true, use CG/PCG with AMG preconditioner.
     *                          If false, use AMG as standalone solver.
     */
    explicit HYPREBoomerAMGSolver(bool useOuterIteration = true);

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

    /**
     * @brief Set the relaxation type for AMG.
     * @param type 0=Jacobi, 1=Gauss-Seidel, 2=SOR, 3=Chebyshev.
     */
    void setRelaxType(int type);

    /**
     * @brief Set the coarsening type for AMG.
     * @param type 0=CLJP, 1=Ruge-Stueben, 2=Falgout, etc.
     */
    void setCoarseningType(int type);

private:
    bool useOuterIteration_;
    int maxIterations_;
    double relativeTolerance_;
    double absoluteTolerance_;
    int printLevel_;
    int relaxType_;
    int coarseningType_;
    mutable int numIterations_;
    mutable double finalResidual_;
};

} // namespace mpfem

#endif
