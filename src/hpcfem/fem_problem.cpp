/**
 * @file fem_problem.cpp
 * @brief Implementation of FemProblem orchestrator
 */

#include "fem_problem.hpp"
#include <stdexcept>

namespace hpcfem
{

#ifdef MFEM_USE_MPI
FemProblem::FemProblem(mfem::ParMesh* pmesh,
                      PhysicsInterface* physics,
                      SolverInterface* solver)
    : pmesh_(pmesh),
      physics_(physics),
      solver_(solver),
      A_(nullptr),
      solutionGf_(nullptr),
      assembled_(false),
      solved_(false)
{
    if (!pmesh_ || !physics_ || !solver_)
    {
        throw std::invalid_argument("FemProblem: null pointer argument");
    }
}
#else
FemProblem::FemProblem(mfem::Mesh* mesh,
                      PhysicsInterface* physics,
                      SolverInterface* solver)
    : mesh_(mesh),
      physics_(physics),
      solver_(solver),
      A_(nullptr),
      solutionGf_(nullptr),
      assembled_(false),
      solved_(false)
{
    if (!mesh_ || !physics_ || !solver_)
    {
        throw std::invalid_argument("FemProblem: null pointer argument");
    }
}
#endif

FemProblem::~FemProblem()
{
    delete A_;
    delete solutionGf_;
}

void FemProblem::assemble()
{
    // Clean up previous assembly
    delete A_;
    A_ = nullptr;
    
    // Initialize vectors
    auto* fespace = physics_->getFiniteElementSpace();
    int size = fespace->GetTrueVSize();
    b_.SetSize(size);
    x_.SetSize(size);
    x_ = 0.0;
    
#ifdef MFEM_USE_MPI
    A_ = new mfem::HypreParMatrix();
#else
    A_ = new mfem::SparseMatrix();
#endif
    
    // Call physics to assemble
    physics_->assemble(*A_, b_, x_, essTdofList_);
    
    assembled_ = true;
    solved_ = false;
}

void FemProblem::solve()
{
    if (!assembled_)
    {
        throw std::runtime_error("FemProblem::solve() called before assemble()");
    }
    
    // Solve the system
    solver_->solve(*A_, b_, x_);
    
    solved_ = true;
}

const mfem::Vector& FemProblem::getSolution() const
{
    if (!solved_)
    {
        throw std::runtime_error("FemProblem::getSolution() called before solve()");
    }
    return x_;
}

#ifdef MFEM_USE_MPI
mfem::ParGridFunction* FemProblem::getSolutionGridFunction()
{
    if (!solved_)
    {
        throw std::runtime_error("FemProblem::getSolutionGridFunction() called before solve()");
    }
    
    if (!solutionGf_)
    {
        solutionGf_ = new mfem::ParGridFunction(physics_->getFiniteElementSpace());
    }
    
    // Recover FEM solution from true DOFs
    solutionGf_->SetFromTrueDofs(x_);
    
    return solutionGf_;
}
#else
mfem::GridFunction* FemProblem::getSolutionGridFunction()
{
    if (!solved_)
    {
        throw std::runtime_error("FemProblem::getSolutionGridFunction() called before solve()");
    }
    
    if (!solutionGf_)
    {
        solutionGf_ = new mfem::GridFunction(physics_->getFiniteElementSpace());
    }
    
    // Recover FEM solution from true DOFs
    solutionGf_->SetFromTrueDofs(x_);
    
    return solutionGf_;
}
#endif

} // namespace hpcfem
