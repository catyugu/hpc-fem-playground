/**
 * @file physics_maxwell_timedomain.cpp
 * @brief Implementation of PhysicsMaxwellTimeDomain class
 */

#include "physics_maxwell_timedomain.hpp"
#include <cmath>
#include <iostream>

namespace hpcfem
{

#ifdef MFEM_USE_MPI

PhysicsMaxwellTimeDomain::PhysicsMaxwellTimeDomain(
    mfem::ParMesh* pmesh,
    int order,
    double (*epsilonFunc)(const mfem::Vector&),
    double (*muInvFunc)(const mfem::Vector&),
    double (*sigmaFunc)(const mfem::Vector&),
    void (*jFunc)(const mfem::Vector&, double, mfem::Vector&),
    const mfem::Array<int>& abcMarkers,
    const mfem::Array<int>& dbcMarkers,
    void (*dEdtFunc)(const mfem::Vector&, double, mfem::Vector&))
    : mfem::TimeDependentOperator(),
      order_(order),
      isLossy_(sigmaFunc != nullptr || abcMarkers.Size() > 0),
      pmesh_(pmesh),
      hCurlFESpace_(nullptr),
      hDivFESpace_(nullptr),
      eField_(nullptr),
      bField_(nullptr),
      dEdtField_(nullptr),
      hDivMassMuInv_(nullptr),
      hCurlLosses_(nullptr),
      weakCurlMuInv_(nullptr),
      curlOp_(nullptr),
      massMuInvMatrix_(nullptr),
      lossMatrix_(nullptr),
      negCurlMatrix_(nullptr),
      weakCurlMatrix_(nullptr),
      eVec_(nullptr),
      bVec_(nullptr),
      rhsVec_(nullptr),
      hCurlFEC_(nullptr),
      hDivFEC_(nullptr),
      materials_(nullptr),
      sources_(nullptr),
      boundaries_(nullptr),
      implicitSolver_(nullptr),
      preconditioner_(nullptr),
      maxTimeStep_(-1.0),
      myRank_(0),
      numProcs_(1)
{
    MPI_Comm_rank(pmesh_->GetComm(), &myRank_);
    MPI_Comm_size(pmesh_->GetComm(), &numProcs_);

    // Create finite element spaces
    hCurlFEC_ = new mfem::ND_FECollection(order_, pmesh_->Dimension());
    hDivFEC_ = new mfem::RT_FECollection(order_, pmesh_->Dimension());

    hCurlFESpace_ = new mfem::ParFiniteElementSpace(pmesh_, hCurlFEC_);
    hDivFESpace_ = new mfem::ParFiniteElementSpace(pmesh_, hDivFEC_);

    // Set operator dimensions
    height = hCurlFESpace_->GlobalTrueVSize();
    width = hDivFESpace_->GlobalTrueVSize();

    // Initialize material properties
    materials_ = new MaxwellMaterialProperties(epsilonFunc, muInvFunc, sigmaFunc);

    // Initialize source terms
    sources_ = new MaxwellSourceTerms(jFunc, hCurlFESpace_);

    // Initialize boundary conditions
    boundaries_ = new MaxwellBoundaryConditions(
        abcMarkers, dbcMarkers, dEdtFunc,
        materials_->getImpedanceInvCoefficient(), hCurlFESpace_);

    // Create grid functions
    eField_ = new mfem::ParGridFunction(hCurlFESpace_);
    bField_ = new mfem::ParGridFunction(hDivFESpace_);
    dEdtField_ = new mfem::ParGridFunction(hCurlFESpace_);

    *eField_ = 0.0;
    *bField_ = 0.0;
    *dEdtField_ = 0.0;

    // Create parallel vectors
    eVec_ = eField_->ParallelProject();
    bVec_ = bField_->ParallelProject();
    rhsVec_ = new mfem::HypreParVector(hCurlFESpace_);

    // Assemble system matrices
    assembleMatrices();

    // Compute maximum time step
    maxTimeStep_ = getMaximumTimeStep();

    if (myRank_ == 0)
    {
        std::cout << "Maxwell solver initialized:" << std::endl;
        std::cout << "  H(curl) DOFs: " << hCurlFESpace_->GlobalTrueVSize() << std::endl;
        std::cout << "  H(div) DOFs: " << hDivFESpace_->GlobalTrueVSize() << std::endl;
        std::cout << "  Max time step: " << maxTimeStep_ << " s" << std::endl;
    }
}

#else

PhysicsMaxwellTimeDomain::PhysicsMaxwellTimeDomain(
    mfem::Mesh* mesh,
    int order,
    double (*epsilonFunc)(const mfem::Vector&),
    double (*muInvFunc)(const mfem::Vector&),
    double (*sigmaFunc)(const mfem::Vector&),
    void (*jFunc)(const mfem::Vector&, double, mfem::Vector&),
    const mfem::Array<int>& abcMarkers,
    const mfem::Array<int>& dbcMarkers,
    void (*dEdtFunc)(const mfem::Vector&, double, mfem::Vector&))
    : mfem::TimeDependentOperator(),
      order_(order),
      isLossy_(sigmaFunc != nullptr || abcMarkers.Size() > 0),
      mesh_(mesh),
      hCurlFESpace_(nullptr),
      hDivFESpace_(nullptr),
      eField_(nullptr),
      bField_(nullptr),
      dEdtField_(nullptr),
      hDivMassMuInv_(nullptr),
      hCurlLosses_(nullptr),
      weakCurlMuInv_(nullptr),
      curlOp_(nullptr),
      massMuInvMatrix_(nullptr),
      lossMatrix_(nullptr),
      negCurlMatrix_(nullptr),
      weakCurlMatrix_(nullptr),
      hCurlFEC_(nullptr),
      hDivFEC_(nullptr),
      materials_(nullptr),
      sources_(nullptr),
      boundaries_(nullptr),
      implicitSolver_(nullptr),
      preconditioner_(nullptr),
      maxTimeStep_(-1.0),
      myRank_(0),
      numProcs_(1)
{
    // Create finite element spaces
    hCurlFEC_ = new mfem::ND_FECollection(order_, mesh_->Dimension());
    hDivFEC_ = new mfem::RT_FECollection(order_, mesh_->Dimension());

    hCurlFESpace_ = new mfem::FiniteElementSpace(mesh_, hCurlFEC_);
    hDivFESpace_ = new mfem::FiniteElementSpace(mesh_, hDivFEC_);

    // Set operator dimensions
    height = hCurlFESpace_->GetTrueVSize();
    width = hDivFESpace_->GetTrueVSize();

    // Initialize material properties
    materials_ = new MaxwellMaterialProperties(epsilonFunc, muInvFunc, sigmaFunc);

    // Initialize source terms
    sources_ = new MaxwellSourceTerms(jFunc, hCurlFESpace_);

    // Initialize boundary conditions
    boundaries_ = new MaxwellBoundaryConditions(
        abcMarkers, dbcMarkers, dEdtFunc,
        materials_->getImpedanceInvCoefficient(), hCurlFESpace_);

    // Create grid functions
    eField_ = new mfem::GridFunction(hCurlFESpace_);
    bField_ = new mfem::GridFunction(hDivFESpace_);
    dEdtField_ = new mfem::GridFunction(hCurlFESpace_);

    *eField_ = 0.0;
    *bField_ = 0.0;
    *dEdtField_ = 0.0;

    // Assemble system matrices
    assembleMatrices();

    // Compute maximum time step
    maxTimeStep_ = getMaximumTimeStep();

    std::cout << "Maxwell solver initialized:" << std::endl;
    std::cout << "  H(curl) DOFs: " << hCurlFESpace_->GetTrueVSize() << std::endl;
    std::cout << "  H(div) DOFs: " << hDivFESpace_->GetTrueVSize() << std::endl;
    std::cout << "  Max time step: " << maxTimeStep_ << " s" << std::endl;
}

#endif

PhysicsMaxwellTimeDomain::~PhysicsMaxwellTimeDomain()
{
    delete implicitSolver_;
    delete preconditioner_;
    delete boundaries_;
    delete sources_;
    delete materials_;

#ifdef MFEM_USE_MPI
    delete rhsVec_;
    delete bVec_;
    delete eVec_;
    delete massMuInvMatrix_;
    delete lossMatrix_;
    delete negCurlMatrix_;
    delete weakCurlMatrix_;
#else
    delete massMuInvMatrix_;
    delete lossMatrix_;
    delete negCurlMatrix_;
    delete weakCurlMatrix_;
#endif

    delete dEdtField_;
    delete bField_;
    delete eField_;
    delete curlOp_;
    delete weakCurlMuInv_;
    delete hCurlLosses_;
    delete hDivMassMuInv_;
    delete hDivFESpace_;
    delete hCurlFESpace_;
    delete hDivFEC_;
    delete hCurlFEC_;
}

void PhysicsMaxwellTimeDomain::assembleMatrices()
{
#ifdef MFEM_USE_MPI
    // H(div) mass matrix with μ⁻¹
    hDivMassMuInv_ = new mfem::ParBilinearForm(hDivFESpace_);
    hDivMassMuInv_->AddDomainIntegrator(
        new mfem::VectorFEMassIntegrator(*materials_->getPermeabilityInvCoefficient()));
    hDivMassMuInv_->Assemble();
    hDivMassMuInv_->Finalize();
    massMuInvMatrix_ = hDivMassMuInv_->ParallelAssemble();

    // Weak curl operator: (μ⁻¹ B, ∇×F)
    weakCurlMuInv_ = new mfem::ParMixedBilinearForm(hDivFESpace_, hCurlFESpace_);
    weakCurlMuInv_->AddDomainIntegrator(
        new mfem::MixedVectorWeakCurlIntegrator(*materials_->getPermeabilityInvCoefficient()));
    weakCurlMuInv_->Assemble();
    weakCurlMuInv_->Finalize();
    weakCurlMatrix_ = weakCurlMuInv_->ParallelAssemble();

    // Curl operator: ∇×E
    curlOp_ = new mfem::ParDiscreteLinearOperator(hCurlFESpace_, hDivFESpace_);
    curlOp_->AddDomainInterpolator(new mfem::CurlInterpolator());
    curlOp_->Assemble();
    curlOp_->Finalize();
    negCurlMatrix_ = curlOp_->ParallelAssemble();
    *negCurlMatrix_ *= -1.0;

    // Loss matrix (if lossy)
    if (isLossy_)
    {
        hCurlLosses_ = new mfem::ParBilinearForm(hCurlFESpace_);
        if (materials_->isLossy())
        {
            hCurlLosses_->AddDomainIntegrator(
                new mfem::VectorFEMassIntegrator(*materials_->getConductivityCoefficient()));
        }
        if (boundaries_->hasAbsorbingBC())
        {
            // ABC contributes to loss matrix
            mfem::ParBilinearForm* abcForm = boundaries_->getAbsorbingForm();
            // TODO: Add ABC form to loss matrix
        }
        hCurlLosses_->Assemble();
        hCurlLosses_->Finalize();
        lossMatrix_ = hCurlLosses_->ParallelAssemble();
    }
#else
    // Serial version
    hDivMassMuInv_ = new mfem::BilinearForm(hDivFESpace_);
    hDivMassMuInv_->AddDomainIntegrator(
        new mfem::VectorFEMassIntegrator(*materials_->getPermeabilityInvCoefficient()));
    hDivMassMuInv_->Assemble();
    hDivMassMuInv_->Finalize();
    massMuInvMatrix_ = &(hDivMassMuInv_->SpMat());

    weakCurlMuInv_ = new mfem::MixedBilinearForm(hDivFESpace_, hCurlFESpace_);
    weakCurlMuInv_->AddDomainIntegrator(
        new mfem::MixedVectorWeakCurlIntegrator(*materials_->getPermeabilityInvCoefficient()));
    weakCurlMuInv_->Assemble();
    weakCurlMuInv_->Finalize();
    weakCurlMatrix_ = &(weakCurlMuInv_->SpMat());

    curlOp_ = new mfem::DiscreteLinearOperator(hCurlFESpace_, hDivFESpace_);
    curlOp_->AddDomainInterpolator(new mfem::CurlInterpolator());
    curlOp_->Assemble();
    curlOp_->Finalize();
    negCurlMatrix_ = &(curlOp_->SpMat());
    *negCurlMatrix_ *= -1.0;

    if (isLossy_)
    {
        hCurlLosses_ = new mfem::BilinearForm(hCurlFESpace_);
        if (materials_->isLossy())
        {
            hCurlLosses_->AddDomainIntegrator(
                new mfem::VectorFEMassIntegrator(*materials_->getConductivityCoefficient()));
        }
        hCurlLosses_->Assemble();
        hCurlLosses_->Finalize();
        lossMatrix_ = &(hCurlLosses_->SpMat());
    }
#endif
}

void PhysicsMaxwellTimeDomain::Mult(const mfem::Vector& B, mfem::Vector& dEdt) const
{
    // Compute dE/dt = (1/ε)[∇×(μ⁻¹B) - σE - J]
    // This is the explicit operator evaluation for lossless case
    
#ifdef MFEM_USE_MPI
    // Step 1: Compute curl of magnetic flux: ∇×(μ⁻¹B)
    // WeakCurlMatrix implements (μ⁻¹B, ∇×F)
    weakCurlMatrix_->Mult(B, *rhsVec_);
    
    // Step 2: Add source term if present
    if (sources_->hasSource())
    {
        mfem::HypreParVector* jDual = sources_->getDualForm()->ParallelAssemble();
        *rhsVec_ -= *jDual;
        delete jDual;
    }
    
    // Step 3: For explicit case, we need to solve M_ε * dE/dt = RHS
    // For now, use a simple Richardson iteration (identity preconditioner)
    // TODO: Use proper mass matrix solver
    
    // For explicit scheme without losses, just return the RHS
    // (assuming mass matrix is identity or using lumped mass)
    dEdt = *rhsVec_;
    
    // Step 4: Apply boundary conditions
    boundaries_->applyDirichletBC(dEdt);
#else
    // Serial version
    weakCurlMatrix_->Mult(B, dEdt);
    
    if (sources_->hasSource())
    {
        mfem::Vector jDual;
        sources_->getDualForm()->Assemble();
        sources_->getDualForm()->ParallelAssemble(jDual);
        dEdt -= jDual;
    }
    
    boundaries_->applyDirichletBC(dEdt);
#endif
}

void PhysicsMaxwellTimeDomain::implicitSolve(double dt, const mfem::Vector& B, 
                                            mfem::Vector& dEdt)
{
    // Implicit solve for lossy case
    // System: (M_ε + dt*M_σ + dt*ABC) * dE/dt = ∇×(μ⁻¹B) - J
    
#ifdef MFEM_USE_MPI
    // Setup solver if not already done
    if (implicitSolver_ == nullptr)
    {
        setupImplicitSolver(dt);
    }
    
    // Compute RHS: ∇×(μ⁻¹B) - J
    weakCurlMatrix_->Mult(B, *rhsVec_);
    
    if (sources_->hasSource())
    {
        mfem::HypreParVector* jDual = sources_->getDualForm()->ParallelAssemble();
        *rhsVec_ -= *jDual;
        delete jDual;
    }
    
    // For now, use explicit method even for lossy case
    // TODO: Assemble system matrix (M_ε + dt*M_σ + dt*ABC)
    // TODO: Solve linear system
    
    // Temporary: just use explicit (will be inaccurate for lossy)
    dEdt = *rhsVec_;
    
    // Apply boundary conditions
    boundaries_->applyDirichletBC(dEdt);
    
    if (myRank_ == 0 && isLossy_)
    {
        std::cout << "Warning: Implicit solve not fully implemented, "
                  << "using explicit approximation" << std::endl;
    }
#else
    // Serial version
    weakCurlMatrix_->Mult(B, dEdt);
    
    if (sources_->hasSource())
    {
        mfem::Vector jDual;
        sources_->getDualForm()->Assemble();
        sources_->getDualForm()->ParallelAssemble(jDual);
        dEdt -= jDual;
    }
    
    boundaries_->applyDirichletBC(dEdt);
#endif
}

double PhysicsMaxwellTimeDomain::getMaximumTimeStep() const
{
    // CFL condition: dt < h / (c * order)
    // where c = 1/sqrt(ε*μ) is speed of light
    
    // Compute minimum mesh size
    double hMin = 1e10;
    
#ifdef MFEM_USE_MPI
    for (int i = 0; i < pmesh_->GetNE(); i++)
    {
        double h = pmesh_->GetElementSize(i);
        if (h < hMin) hMin = h;
    }
    double hMinGlobal;
    MPI_Allreduce(&hMin, &hMinGlobal, 1, MPI_DOUBLE, MPI_MIN, pmesh_->GetComm());
    hMin = hMinGlobal;
#else
    for (int i = 0; i < mesh_->GetNE(); i++)
    {
        double h = mesh_->GetElementSize(i);
        if (h < hMin) hMin = h;
    }
#endif

    // Speed of light in free space
    const double c0 = 2.99792458e8;  // m/s
    
    // Safety factor
    const double safetyFactor = 0.5;
    
    return safetyFactor * hMin / (c0 * order_);
}

double PhysicsMaxwellTimeDomain::getEnergy() const
{
    // Compute electromagnetic energy: E = 0.5 * (E^T M_ε E + B^T M_μ B)
    
    double energyE = 0.0;
    double energyB = 0.0;
    
#ifdef MFEM_USE_MPI
    // Electric field energy: 0.5 * E^T M_ε E
    // For simplicity, approximate with dot product (assumes identity mass matrix)
    // TODO: Use actual permittivity mass matrix
    energyE = 0.5 * mfem::InnerProduct(*eVec_, *eVec_);
    
    // Magnetic field energy: 0.5 * B^T M_μ B
    mfem::HypreParVector temp(hDivFESpace_);
    massMuInvMatrix_->Mult(*bVec_, temp);
    energyB = 0.5 * mfem::InnerProduct(*bVec_, temp);
    
    // Sum across processors
    double localEnergy = energyE + energyB;
    double globalEnergy = 0.0;
    MPI_Allreduce(&localEnergy, &globalEnergy, 1, MPI_DOUBLE, MPI_SUM, 
                  pmesh_->GetComm());
    return globalEnergy;
#else
    // Serial version
    energyE = 0.5 * (*eField_) * (*eField_);
    
    mfem::Vector temp(hDivFESpace_->GetTrueVSize());
    massMuInvMatrix_->Mult(*bField_, temp);
    energyB = 0.5 * (*bField_) * temp;
    
    return energyE + energyB;
#endif
}

void PhysicsMaxwellTimeDomain::setInitialEField(mfem::VectorCoefficient& eCoef)
{
    eField_->ProjectCoefficient(eCoef);
#ifdef MFEM_USE_MPI
    eField_->ParallelProject(*eVec_);
#endif
}

void PhysicsMaxwellTimeDomain::setInitialBField(mfem::VectorCoefficient& bCoef)
{
    bField_->ProjectCoefficient(bCoef);
#ifdef MFEM_USE_MPI
    bField_->ParallelProject(*bVec_);
#endif
}

void PhysicsMaxwellTimeDomain::syncGridFunctions()
{
#ifdef MFEM_USE_MPI
    eField_->Distribute(eVec_);
    bField_->Distribute(bVec_);
#endif
}

void PhysicsMaxwellTimeDomain::setupImplicitSolver(double dt)
{
    // Setup linear solver for implicit time integration
    // System: (M_ε + dt*M_σ + dt*ABC) * dE/dt = RHS
    
#ifdef MFEM_USE_MPI
    if (implicitSolver_ != nullptr)
    {
        delete implicitSolver_;
        delete preconditioner_;
    }
    
    // Create CG solver
    implicitSolver_ = new mfem::CGSolver(pmesh_->GetComm());
    implicitSolver_->SetRelTol(1e-8);
    implicitSolver_->SetMaxIter(1000);
    implicitSolver_->SetPrintLevel(0);
    
    // Create diagonal preconditioner (Jacobi)
    preconditioner_ = new mfem::HypreSmoother();
    preconditioner_->SetType(mfem::HypreSmoother::Jacobi);
    
    implicitSolver_->SetPreconditioner(*preconditioner_);
    
    if (myRank_ == 0)
    {
        std::cout << "Implicit solver configured for dt = " << dt << std::endl;
    }
#else
    // Serial version - use simple CG
    if (implicitSolver_ != nullptr)
    {
        delete implicitSolver_;
    }
    
    implicitSolver_ = new mfem::CGSolver();
    implicitSolver_->SetRelTol(1e-8);
    implicitSolver_->SetMaxIter(1000);
    implicitSolver_->SetPrintLevel(0);
#endif
}

} // namespace hpcfem
