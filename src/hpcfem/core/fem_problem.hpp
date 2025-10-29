/**
 * @file core/fem_problem.hpp
 * @brief Top-level FEM problem orchestrator (moved to core/)
 */

#ifndef HPCFEM_CORE_FEM_PROBLEM_HPP
#define HPCFEM_CORE_FEM_PROBLEM_HPP

#include "hpcfem/core/physics_interface.hpp"
#include "hpcfem/core/solver_interface.hpp"
#include "mfem.hpp"

namespace hpcfem
{

class FemProblem
{
public:
#ifdef MFEM_USE_MPI
    FemProblem(mfem::ParMesh* pmesh,
               PhysicsInterface* physics,
               SolverInterface* solver);
#else
    FemProblem(mfem::Mesh* mesh,
               PhysicsInterface* physics,
               SolverInterface* solver);
#endif

    ~FemProblem();

    void assemble();
    void solve();
    const mfem::Vector& getSolution() const;

#ifdef MFEM_USE_MPI
    mfem::ParGridFunction* getSolutionGridFunction();
#else
    mfem::GridFunction* getSolutionGridFunction();
#endif

private:
#ifdef MFEM_USE_MPI
    mfem::ParMesh* pmesh_;
    mfem::HypreParMatrix* A_;
    mfem::ParGridFunction* solutionGf_;
#else
    mfem::Mesh* mesh_;
    mfem::SparseMatrix* A_;
    mfem::GridFunction* solutionGf_;
#endif

    PhysicsInterface* physics_;
    SolverInterface* solver_;
    
    mfem::Vector b_;
    mfem::Vector x_;
    mfem::Array<int> essTdofList_;
    
    bool assembled_;
    bool solved_;
};

} // namespace hpcfem

#endif // HPCFEM_CORE_FEM_PROBLEM_HPP
