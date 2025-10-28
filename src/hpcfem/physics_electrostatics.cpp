/**
 * @file physics_electrostatics.cpp
 * @brief Implementation of electrostatics physics
 */

#include "physics_electrostatics.hpp"

namespace hpcfem
{

#ifdef MFEM_USE_MPI
ElectrostaticsPhysics::ElectrostaticsPhysics(mfem::ParMesh* pmesh,
                                           int polynomialOrder,
                                           mfem::Coefficient* conductivity,
                                           mfem::Coefficient* sourceCoeff,
                                           mfem::Coefficient* boundaryCoeff)
    : conductivity_(conductivity),
      sourceCoeff_(sourceCoeff),
      boundaryCoeff_(boundaryCoeff)
{
    int dim = pmesh->Dimension();
    fec_ = new mfem::H1_FECollection(polynomialOrder, dim);
    fespace_ = new mfem::ParFiniteElementSpace(pmesh, fec_);
    
    bilinearForm_ = new mfem::ParBilinearForm(fespace_);
    bilinearForm_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*conductivity_));
    
    linearForm_ = new mfem::ParLinearForm(fespace_);
    linearForm_->AddDomainIntegrator(new mfem::DomainLFIntegrator(*sourceCoeff_));
}
#else
ElectrostaticsPhysics::ElectrostaticsPhysics(mfem::Mesh* mesh,
                                           int polynomialOrder,
                                           mfem::Coefficient* conductivity,
                                           mfem::Coefficient* sourceCoeff,
                                           mfem::Coefficient* boundaryCoeff)
    : conductivity_(conductivity),
      sourceCoeff_(sourceCoeff),
      boundaryCoeff_(boundaryCoeff)
{
    int dim = mesh->Dimension();
    fec_ = new mfem::H1_FECollection(polynomialOrder, dim);
    fespace_ = new mfem::FiniteElementSpace(mesh, fec_);
    
    bilinearForm_ = new mfem::BilinearForm(fespace_);
    bilinearForm_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*conductivity_));
    
    linearForm_ = new mfem::LinearForm(fespace_);
    linearForm_->AddDomainIntegrator(new mfem::DomainLFIntegrator(*sourceCoeff_));
}
#endif

ElectrostaticsPhysics::~ElectrostaticsPhysics()
{
    delete linearForm_;
    delete bilinearForm_;
    delete fespace_;
    delete fec_;
}

#ifdef MFEM_USE_MPI
void ElectrostaticsPhysics::assemble(mfem::HypreParMatrix& A,
                                    mfem::Vector& b,
                                    mfem::Vector& x,
                                    mfem::Array<int>& essTdofList)
{
    bilinearForm_->Assemble();
    linearForm_->Assemble();
    
    // Apply boundary conditions
    mfem::ParGridFunction boundaryFunc(fespace_);
    boundaryFunc.ProjectCoefficient(*boundaryCoeff_);
    
    // Get essential boundary DOFs (all boundaries)
    mfem::Array<int> essBdr(fespace_->GetParMesh()->bdr_attributes.Max());
    essBdr = 1;
    fespace_->GetEssentialTrueDofs(essBdr, essTdofList);
    
    // Form linear system
    bilinearForm_->FormLinearSystem(essTdofList, boundaryFunc, *linearForm_, A, x, b);
}

mfem::ParFiniteElementSpace* ElectrostaticsPhysics::getFiniteElementSpace()
{
    return fespace_;
}
#else
void ElectrostaticsPhysics::assemble(mfem::SparseMatrix& A,
                                    mfem::Vector& b,
                                    mfem::Vector& x,
                                    mfem::Array<int>& essTdofList)
{
    bilinearForm_->Assemble();
    linearForm_->Assemble();
    
    // Apply boundary conditions
    mfem::GridFunction boundaryFunc(fespace_);
    boundaryFunc.ProjectCoefficient(*boundaryCoeff_);
    
    // Get essential boundary DOFs (all boundaries)
    mfem::Array<int> essBdr(fespace_->GetMesh()->bdr_attributes.Max());
    essBdr = 1;
    fespace_->GetEssentialTrueDofs(essBdr, essTdofList);
    
    // Form linear system
    bilinearForm_->FormLinearSystem(essTdofList, boundaryFunc, *linearForm_, A, x, b);
}

mfem::FiniteElementSpace* ElectrostaticsPhysics::getFiniteElementSpace()
{
    return fespace_;
}
#endif

} // namespace hpcfem
