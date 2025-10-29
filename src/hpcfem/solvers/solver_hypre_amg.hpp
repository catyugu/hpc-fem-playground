/**
 * @file solvers/solver_hypre_amg.hpp
 * @brief Concrete HYPRE BoomerAMG solver implementation
 */

#ifndef HPCFEM_SOLVER_HYPRE_AMG_HPP
#define HPCFEM_SOLVER_HYPRE_AMG_HPP

#include "hpcfem/core/solver_interface.hpp"
#include "mfem.hpp"

namespace hpcfem
{

class HypreAmgSolver : public SolverInterface
{
public:
    HypreAmgSolver(double relTol = 1.0e-12, int maxIter = 2000, int printLevel = 0);
    ~HypreAmgSolver() override = default;

#ifdef MFEM_USE_MPI
    void solve(const mfem::HypreParMatrix& A, const mfem::Vector& b, mfem::Vector& x) override;
#else
    void solve(const mfem::SparseMatrix& A, const mfem::Vector& b, mfem::Vector& x) override;
#endif

    int getNumIterations() const { return numIterations_; }
    double getFinalNorm() const { return finalNorm_; }

private:
    double relTol_;
    int maxIter_;
    int printLevel_;
    int numIterations_;
    double finalNorm_;
};

} // namespace hpcfem

#endif // HPCFEM_SOLVER_HYPRE_AMG_HPP
