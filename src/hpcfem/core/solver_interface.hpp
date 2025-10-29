/**
 * @file core/solver_interface.hpp
 * @brief Abstract base class for all linear solvers in the hpcfem library
 */

#ifndef HPCFEM_CORE_SOLVER_INTERFACE_HPP
#define HPCFEM_CORE_SOLVER_INTERFACE_HPP

#include "mfem.hpp"

namespace hpcfem
{

class SolverInterface
{
public:
    virtual ~SolverInterface() = default;

#ifdef MFEM_USE_MPI
    virtual void solve(const mfem::HypreParMatrix& A, const mfem::Vector& b, mfem::Vector& x) = 0;
#else
    virtual void solve(const mfem::SparseMatrix& A, const mfem::Vector& b, mfem::Vector& x) = 0;
#endif
};

} // namespace hpcfem

#endif // HPCFEM_CORE_SOLVER_INTERFACE_HPP
