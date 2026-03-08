#ifndef MPFEM_LINEAR_SYSTEM_SOLVER_HPP
#define MPFEM_LINEAR_SYSTEM_SOLVER_HPP

#include "mfem.hpp"

#include <string>

namespace mpfem {

/**
 * @brief Generic linear solve wrapper for SPD systems assembled by MFEM.
 */
class LinearSystemSolver {
public:
    /**
     * @brief Solves A x = b using CG with Gauss-Seidel preconditioner.
     * @param matrix System matrix.
     * @param rhs Right-hand side.
     * @param solution Initial guess and output solution.
     * @param maxIterations Maximum number of iterations.
     * @param relativeTolerance Relative residual tolerance.
     * @param errorMessage Error details on failure.
     * @return True when solver succeeds.
     */
    static bool solve(mfem::SparseMatrix &matrix,
                      mfem::Vector &rhs,
                      mfem::Vector &solution,
                      int maxIterations,
                      double relativeTolerance,
                      std::string &errorMessage);
};

} // namespace mpfem

#endif
