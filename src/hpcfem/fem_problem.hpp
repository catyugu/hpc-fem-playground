/**
 * @file fem_problem.hpp
 * @brief Top-level FEM problem orchestrator
 * 
 * This class coordinates mesh, physics, and solver to provide
 * a unified interface for solving FEM problems.
 */

#ifndef HPCFEM_FEM_PROBLEM_HPP
#define HPCFEM_FEM_PROBLEM_HPP

#include "physics_interface.hpp"
#include "solver_interface.hpp"
#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class FemProblem
 * @brief Orchestrates the complete FEM problem workflow
 * 
 * This class owns or references all components needed to solve a FEM problem:
 * - Mesh (referenced, not owned)
 * - Physics (referenced, not owned)
 * - Solver (referenced, not owned)
 * 
 * It provides high-level methods:
 * - assemble(): Build the linear system
 * - solve(): Solve the linear system
 * - getSolution(): Retrieve the solution vector
 */
class FemProblem
{
public:
#ifdef MFEM_USE_MPI
    /**
     * @brief Constructor for parallel execution
     * @param pmesh Parallel mesh (not owned)
     * @param physics Physics module (not owned)
     * @param solver Linear solver (not owned)
     */
    FemProblem(mfem::ParMesh* pmesh,
              PhysicsInterface* physics,
              SolverInterface* solver);
#else
    /**
     * @brief Constructor for serial execution
     * @param mesh Mesh (not owned)
     * @param physics Physics module (not owned)
     * @param solver Linear solver (not owned)
     */
    FemProblem(mfem::Mesh* mesh,
              PhysicsInterface* physics,
              SolverInterface* solver);
#endif

    /**
     * @brief Destructor
     */
    ~FemProblem();

    /**
     * @brief Assemble the linear system
     * 
     * Calls physics->assemble() to build A, b, and apply boundary conditions.
     */
    void assemble();

    /**
     * @brief Solve the linear system
     * 
     * Calls solver->solve() to compute the solution.
     */
    void solve();

    /**
     * @brief Get the solution vector
     * @return Reference to the solution vector
     */
    const mfem::Vector& getSolution() const;

    /**
     * @brief Get the solution as a grid function
     * @return Pointer to grid function (ownership retained)
     */
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

#endif // HPCFEM_FEM_PROBLEM_HPP
