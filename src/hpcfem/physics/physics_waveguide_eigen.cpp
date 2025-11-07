/**
 * @file physics/physics_waveguide_eigen.cpp
 * @brief Implementation of 3D Maxwell eigenvalue problem solver
 * 
 * Based on MFEM example ex13p.cpp
 */

#include "physics_waveguide_eigen.hpp"
#include <iostream>
#include <stdexcept>

namespace hpcfem
{

#ifdef MFEM_USE_MPI
PhysicsWaveguideEigen::PhysicsWaveguideEigen(mfem::ParMesh* pmesh, int order)
    : pmesh_(pmesh),
      curlCurlForm_(nullptr),
      massForm_(nullptr),
      fespace_(nullptr),
      feCollection_(nullptr),
      polynomialOrder_(order),
      isSolved_(false)
{
    if (!pmesh_)
    {
        throw std::runtime_error("PhysicsWaveguideEigen: mesh pointer is null");
    }
    
    int dim = pmesh_->Dimension();
    if (dim != 3)
    {
        throw std::runtime_error("PhysicsWaveguideEigen: only 3D meshes are supported");
    }
    
    // Create Nedelec finite element collection (edge elements for E-field)
    feCollection_ = new mfem::ND_FECollection(polynomialOrder_, dim);
    fespace_ = new mfem::ParFiniteElementSpace(pmesh_, feCollection_);
}
#else
PhysicsWaveguideEigen::PhysicsWaveguideEigen(mfem::Mesh* mesh, int order)
    : mesh_(mesh),
      curlCurlForm_(nullptr),
      massForm_(nullptr),
      fespace_(nullptr),
      feCollection_(nullptr),
      polynomialOrder_(order),
      isSolved_(false)
{
    if (!mesh_)
    {
        throw std::runtime_error("PhysicsWaveguideEigen: mesh pointer is null");
    }
    
    int dim = mesh_->Dimension();
    if (dim != 3)
    {
        throw std::runtime_error("PhysicsWaveguideEigen: only 3D meshes are supported");
    }
    
    // Create Nedelec finite element collection (edge elements for E-field)
    feCollection_ = new mfem::ND_FECollection(polynomialOrder_, dim);
    fespace_ = new mfem::FiniteElementSpace(mesh_, feCollection_);
}
#endif

PhysicsWaveguideEigen::~PhysicsWaveguideEigen()
{
    delete curlCurlForm_;
    delete massForm_;
    delete fespace_;
    delete feCollection_;
}

void PhysicsWaveguideEigen::assembleCurlCurlMatrix()
{
    // Assemble the curl-curl bilinear form: (∇×E, ∇×F)
    // This represents the stiffness matrix [A] in the eigenvalue problem
#ifdef MFEM_USE_MPI
    curlCurlForm_ = new mfem::ParBilinearForm(fespace_);
#else
    curlCurlForm_ = new mfem::BilinearForm(fespace_);
#endif
    
    mfem::ConstantCoefficient oneCoeff(1.0);
    curlCurlForm_->AddDomainIntegrator(new mfem::CurlCurlIntegrator(oneCoeff));
    
    // For closed cavities/periodic meshes without boundaries, add mass term
    // to handle the null space (as in ex13p.cpp)
#ifdef MFEM_USE_MPI
    if (pmesh_->bdr_attributes.Size() == 0)
    {
        curlCurlForm_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(oneCoeff));
    }
#else
    if (mesh_->bdr_attributes.Size() == 0)
    {
        curlCurlForm_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(oneCoeff));
    }
#endif
    
    curlCurlForm_->Assemble();
    
    // Apply PEC boundary conditions: E × n = 0
    // Mark all boundary attributes as essential (Dirichlet)
#ifdef MFEM_USE_MPI
    if (pmesh_->bdr_attributes.Size())
    {
        mfem::Array<int> essBdr(pmesh_->bdr_attributes.Max());
        essBdr = 1;
        curlCurlForm_->EliminateEssentialBCDiag(essBdr, 1.0);
    }
#else
    if (mesh_->bdr_attributes.Size())
    {
        mfem::Array<int> essBdr(mesh_->bdr_attributes.Max());
        essBdr = 1;
        curlCurlForm_->EliminateEssentialBCDiag(essBdr, 1.0);
    }
#endif
    
    curlCurlForm_->Finalize();
}

void PhysicsWaveguideEigen::assembleMassMatrix()
{
    // Assemble the mass matrix bilinear form: (E, F)
    // This represents the mass matrix [M] in the eigenvalue problem
#ifdef MFEM_USE_MPI
    massForm_ = new mfem::ParBilinearForm(fespace_);
#else
    massForm_ = new mfem::BilinearForm(fespace_);
#endif
    
    mfem::ConstantCoefficient oneCoeff(1.0);
    massForm_->AddDomainIntegrator(new mfem::VectorFEMassIntegrator(oneCoeff));
    massForm_->Assemble();
    
    // Shift the eigenvalue corresponding to eliminated DOFs to large value
    // This pushes Dirichlet eigenvalues out of the computational range
#ifdef MFEM_USE_MPI
    if (pmesh_->bdr_attributes.Size())
    {
        mfem::Array<int> essBdr(pmesh_->bdr_attributes.Max());
        essBdr = 1;
        massForm_->EliminateEssentialBCDiag(essBdr, std::numeric_limits<double>::min());
    }
#else
    if (mesh_->bdr_attributes.Size())
    {
        mfem::Array<int> essBdr(mesh_->bdr_attributes.Max());
        essBdr = 1;
        massForm_->EliminateEssentialBCDiag(essBdr, std::numeric_limits<double>::min());
    }
#endif
    
    massForm_->Finalize();
}

std::vector<double> PhysicsWaveguideEigen::solveEigenvalues(int numModes)
{
    if (numModes <= 0)
    {
        throw std::runtime_error("PhysicsWaveguideEigen: numModes must be positive");
    }

    // Clear previous results
    eigenvalues_.clear();
    eigenvectors_.clear();

    // Assemble matrices
    assembleCurlCurlMatrix();
    assembleMassMatrix();

#ifdef MFEM_USE_MPI
    // Get the parallel sparse matrix operators
    mfem::HypreParMatrix* A = curlCurlForm_->ParallelAssemble();
    mfem::HypreParMatrix* M = massForm_->ParallelAssemble();
    
    // Setup HypreAMS preconditioner (Auxiliary-space Maxwell Solver)
    // This is the correct preconditioner for Maxwell eigenvalue problems
    mfem::HypreAMS* ams = new mfem::HypreAMS(*A, fespace_);
    ams->SetPrintLevel(0);
    ams->SetSingularProblem(); // Handle null space of curl-curl operator
    
    // Use HypreAME eigensolver (as in ex13p.cpp)
    // This uses LOBPCG with AMS preconditioner internally
    mfem::HypreAME* ame = new mfem::HypreAME(MPI_COMM_WORLD);
    ame->SetNumModes(numModes);
    ame->SetPreconditioner(*ams);
    ame->SetMaxIter(100);
    ame->SetTol(1e-8);
    ame->SetPrintLevel(0); // Set to 1 for debugging
    ame->SetMassMatrix(*M);
    ame->SetOperator(*A);
    
    // Solve the eigenvalue problem
    ame->Solve();
    
    // Extract eigenvalues and eigenvectors
    mfem::Array<double> evals;
    ame->GetEigenvalues(evals);
    
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0)
    {
        std::cout << "Computed " << evals.Size() << " eigenvalues:" << std::endl;
        for (int i = 0; i < evals.Size(); ++i)
        {
            std::cout << "  Mode " << i << ": λ = " << evals[i] << std::endl;
        }
    }
    
    // Store eigenvalues and eigenvectors (store local DOFs vectors)
    for (int i = 0; i < evals.Size(); ++i)
    {
        eigenvalues_.push_back(evals[i]);

        // Convert HypreParVector (true DOFs) to a ParGridFunction which
        // distributes the true DOFs to local (per-process) DOFs, then
        // copy the local DOF vector into an mfem::Vector and store it.
        mfem::ParGridFunction x(fespace_);
        x = ame->GetEigenvector(i); // Distribute true DOFs into local DOFs

        mfem::Vector localVec; // will contain local DOFs (size = GetVSize())
        localVec = x; // GridFunction is-a Vector; this copies local data

        eigenvectors_.push_back(localVec);
    }
    
    delete ame;
    delete ams;
    delete A;
    delete M;
#else
    // Serial version not implemented - use MPI version
    throw std::runtime_error("PhysicsWaveguideEigen: Serial eigenvalue solver "
                           "not implemented. Please use MPI version.");
#endif

    isSolved_ = true;
    return eigenvalues_;
}

std::vector<mfem::Vector> PhysicsWaveguideEigen::getEigenvectors() const
{
    if (!isSolved_)
    {
        throw std::runtime_error("PhysicsWaveguideEigen: must call solveEigenvalues() first");
    }
    return eigenvectors_;
}

#ifdef MFEM_USE_MPI
mfem::ParGridFunction* PhysicsWaveguideEigen::getModeShape(int modeIndex)
{
    if (!isSolved_)
    {
        throw std::runtime_error("PhysicsWaveguideEigen: must call solveEigenvalues() first");
    }
    
    if (modeIndex < 0 || modeIndex >= static_cast<int>(eigenvectors_.size()))
    {
        throw std::runtime_error("PhysicsWaveguideEigen: mode index out of range");
    }
    
    mfem::ParGridFunction* gf = new mfem::ParGridFunction(fespace_);
    *gf = eigenvectors_[modeIndex];
    return gf;
}
#else
mfem::GridFunction* PhysicsWaveguideEigen::getModeShape(int modeIndex)
{
    if (!isSolved_)
    {
        throw std::runtime_error("PhysicsWaveguideEigen: must call solveEigenvalues() first");
    }
    
    if (modeIndex < 0 || modeIndex >= static_cast<int>(eigenvectors_.size()))
    {
        throw std::runtime_error("PhysicsWaveguideEigen: mode index out of range");
    }
    
    mfem::GridFunction* gf = new mfem::GridFunction(fespace_);
    *gf = eigenvectors_[modeIndex];
    return gf;
}
#endif

} // namespace hpcfem
