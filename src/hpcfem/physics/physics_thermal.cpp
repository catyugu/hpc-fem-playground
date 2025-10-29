/**
 * @file physics/physics_thermal.cpp
 * @brief Implementation of thermal diffusion physics
 */

#include "hpcfem/physics/physics_thermal.hpp"

namespace hpcfem
{

#ifdef MFEM_USE_MPI
ThermalPhysics::ThermalPhysics(mfem::ParMesh* pmesh,
                               int polynomialOrder,
                               mfem::Coefficient* thermalConductivity,
                               mfem::Coefficient* heatSource,
                               mfem::Coefficient* boundaryCoeff)
    : thermalConductivity_(thermalConductivity),
      heatSource_(heatSource),
      boundaryCoeff_(boundaryCoeff)
{
    int dim = pmesh->Dimension();
    fec_ = new mfem::H1_FECollection(polynomialOrder, dim);
    fespace_ = new mfem::ParFiniteElementSpace(pmesh, fec_);
    
    bilinearForm_ = new mfem::ParBilinearForm(fespace_);
    bilinearForm_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*thermalConductivity_));
    
    linearForm_ = new mfem::ParLinearForm(fespace_);
    linearForm_->AddDomainIntegrator(new mfem::DomainLFIntegrator(*heatSource_));
}
#else
ThermalPhysics::ThermalPhysics(mfem::Mesh* mesh,
                               int polynomialOrder,
                               mfem::Coefficient* thermalConductivity,
                               mfem::Coefficient* heatSource,
                               mfem::Coefficient* boundaryCoeff)
    : thermalConductivity_(thermalConductivity),
      heatSource_(heatSource),
      boundaryCoeff_(boundaryCoeff)
{
    int dim = mesh->Dimension();
    fec_ = new mfem::H1_FECollection(polynomialOrder, dim);
    fespace_ = new mfem::FiniteElementSpace(mesh, fec_);
    
    bilinearForm_ = new mfem::BilinearForm(fespace_);
    bilinearForm_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*thermalConductivity_));
    
    linearForm_ = new mfem::LinearForm(fespace_);
    linearForm_->AddDomainIntegrator(new mfem::DomainLFIntegrator(*heatSource_));
}
#endif

ThermalPhysics::~ThermalPhysics()
{
    delete linearForm_;
    delete bilinearForm_;
    delete fespace_;
    delete fec_;
}

#ifdef MFEM_USE_MPI
void ThermalPhysics::assemble(mfem::HypreParMatrix& A,
                              mfem::Vector& b,
                              mfem::Vector& x,
                              mfem::Array<int>& essTdofList)
{
    bilinearForm_->Assemble();
    linearForm_->Assemble();
    
    mfem::ParGridFunction boundaryFunc(fespace_);
    boundaryFunc.ProjectCoefficient(*boundaryCoeff_);
    
    mfem::Array<int> essBdr(fespace_->GetParMesh()->bdr_attributes.Max());
    essBdr = 1;
    fespace_->GetEssentialTrueDofs(essBdr, essTdofList);
    
    bilinearForm_->FormLinearSystem(essTdofList, boundaryFunc, *linearForm_, A, x, b);
}

mfem::ParFiniteElementSpace* ThermalPhysics::getFiniteElementSpace()
{
    return fespace_;
}
#else
void ThermalPhysics::assemble(mfem::SparseMatrix& A,
                              mfem::Vector& b,
                              mfem::Vector& x,
                              mfem::Array<int>& essTdofList)
{
    bilinearForm_->Assemble();
    linearForm_->Assemble();
    
    mfem::GridFunction boundaryFunc(fespace_);
    boundaryFunc.ProjectCoefficient(*boundaryCoeff_);
    
    mfem::Array<int> essBdr(fespace_->GetMesh()->bdr_attributes.Max());
    essBdr = 1;
    fespace_->GetEssentialTrueDofs(essBdr, essTdofList);
    
    bilinearForm_->FormLinearSystem(essTdofList, boundaryFunc, *linearForm_, A, x, b);
}

mfem::FiniteElementSpace* ThermalPhysics::getFiniteElementSpace()
{
    return fespace_;
}
#endif

} // namespace hpcfem
