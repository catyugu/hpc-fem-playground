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
      hCurlMassEpsilon_(nullptr),
      hDivMassMuInv_(nullptr),
      hCurlLosses_(nullptr),
      weakCurlMuInv_(nullptr),
      curlOp_(nullptr),
      massEpsilonMatrix_(nullptr),
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
      hCurlMassEpsilon_(nullptr),
      hDivMassMuInv_(nullptr),
      hCurlLosses_(nullptr),
      weakCurlMuInv_(nullptr),
      curlOp_(nullptr),
      massEpsilonMatrix_(nullptr),
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
    delete massEpsilonMatrix_;
    delete massMuInvMatrix_;
    delete lossMatrix_;
    delete negCurlMatrix_;
    delete weakCurlMatrix_;
#else
    // Serial: some matrices created by Add() need to be deleted
    if (lossMatrix_ != nullptr && lossMatrix_ != &(hCurlLosses_->SpMat()))
    {
        delete lossMatrix_;  // Created by Add() for ABC+conductivity
    }
    // Other matrices are owned by forms, don't delete
#endif

    delete dEdtField_;
    delete bField_;
    delete eField_;
    delete curlOp_;
    delete weakCurlMuInv_;
    delete hCurlLosses_;
    delete hCurlMassEpsilon_;
    delete hDivMassMuInv_;
    delete hDivFESpace_;
    delete hCurlFESpace_;
    delete hDivFEC_;
    delete hCurlFEC_;
}

void PhysicsMaxwellTimeDomain::assembleMatrices()
{
#ifdef MFEM_USE_MPI
    // H(curl) mass matrix with ε
    hCurlMassEpsilon_ = new mfem::ParBilinearForm(hCurlFESpace_);
    hCurlMassEpsilon_->AddDomainIntegrator(
        new mfem::VectorFEMassIntegrator(*materials_->getPermittivityCoefficient()));
    hCurlMassEpsilon_->Assemble();
    hCurlMassEpsilon_->Finalize();
    massEpsilonMatrix_ = hCurlMassEpsilon_->ParallelAssemble();

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

    // Loss matrix (if lossy or has ABC)
    if (isLossy_)
    {
        hCurlLosses_ = new mfem::ParBilinearForm(hCurlFESpace_);
        if (materials_->isLossy())
        {
            hCurlLosses_->AddDomainIntegrator(
                new mfem::VectorFEMassIntegrator(*materials_->getConductivityCoefficient()));
        }
        hCurlLosses_->Assemble();
        hCurlLosses_->Finalize();
        
        mfem::HypreParMatrix* conductivityMatrix = hCurlLosses_->ParallelAssemble();
        
        // Add ABC contribution if present
        if (boundaries_->hasAbsorbingBC())
        {
            // ABC form is already assembled in boundaries_
            mfem::HypreParMatrix* abcMatrix = boundaries_->getAbsorbingForm()->ParallelAssemble();
            
            // Combine conductivity and ABC: M_loss = M_σ + M_ABC
            lossMatrix_ = mfem::Add(1.0, *conductivityMatrix, 1.0, *abcMatrix);
            
            delete conductivityMatrix;
            delete abcMatrix;
        }
        else
        {
            lossMatrix_ = conductivityMatrix;
        }
    }
#else
    // Serial version
    // H(curl) mass matrix with ε
    hCurlMassEpsilon_ = new mfem::BilinearForm(hCurlFESpace_);
    hCurlMassEpsilon_->AddDomainIntegrator(
        new mfem::VectorFEMassIntegrator(*materials_->getPermittivityCoefficient()));
    hCurlMassEpsilon_->Assemble();
    hCurlMassEpsilon_->Finalize();
    massEpsilonMatrix_ = &(hCurlMassEpsilon_->SpMat());

    // H(div) mass matrix with μ⁻¹
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
        
        mfem::SparseMatrix* conductivityMatrix = &(hCurlLosses_->SpMat());
        
        // Add ABC contribution if present
        if (boundaries_->hasAbsorbingBC())
        {
            mfem::SparseMatrix* abcMatrix = &(boundaries_->getAbsorbingForm()->SpMat());
            
            // Combine conductivity and ABC
            lossMatrix_ = mfem::Add(1.0, *conductivityMatrix, 1.0, *abcMatrix);
        }
        else
        {
            lossMatrix_ = conductivityMatrix;
        }
    }
#endif
}

void PhysicsMaxwellTimeDomain::Mult(const mfem::Vector& B, mfem::Vector& dEdt) const
{
    // Compute dE/dt = M_ε⁻¹ * [∇×(μ⁻¹B) - σE - J]
    // Step 1: Compute RHS = ∇×(μ⁻¹B) - J
    
#ifdef MFEM_USE_MPI
    // Compute curl of magnetic flux: ∇×(μ⁻¹B)
    weakCurlMatrix_->Mult(B, *rhsVec_);
    
    // Subtract source term if present
    if (sources_->hasSource())
    {
        mfem::HypreParVector* jDual = sources_->getDualForm()->ParallelAssemble();
        *rhsVec_ -= *jDual;
        delete jDual;
    }
    
    // Step 2: Solve M_ε * dE/dt = RHS using CG
    // Setup solver if not already done
    if (implicitSolver_ == nullptr)
    {
        // Use const_cast since we're setting up solver lazily
        PhysicsMaxwellTimeDomain* nonConstThis = 
            const_cast<PhysicsMaxwellTimeDomain*>(this);
        
        nonConstThis->implicitSolver_ = new mfem::CGSolver(pmesh_->GetComm());
        nonConstThis->implicitSolver_->SetRelTol(1e-12);
        nonConstThis->implicitSolver_->SetMaxIter(1000);
        nonConstThis->implicitSolver_->SetPrintLevel(0);
        
        nonConstThis->preconditioner_ = new mfem::HypreSmoother();
        nonConstThis->preconditioner_->SetType(mfem::HypreSmoother::Jacobi);
        nonConstThis->preconditioner_->SetOperator(*massEpsilonMatrix_);
        nonConstThis->implicitSolver_->SetPreconditioner(*nonConstThis->preconditioner_);
        nonConstThis->implicitSolver_->SetOperator(*massEpsilonMatrix_);
    }
    
    // Solve for dE/dt
    dEdt = 0.0;
    implicitSolver_->Mult(*rhsVec_, dEdt);
    
    // Apply boundary conditions
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
    
    // Solve M_ε * y = dEdt for y
    if (implicitSolver_ == nullptr)
    {
        PhysicsMaxwellTimeDomain* nonConstThis = 
            const_cast<PhysicsMaxwellTimeDomain*>(this);
        
        nonConstThis->implicitSolver_ = new mfem::CGSolver();
        nonConstThis->implicitSolver_->SetOperator(*massEpsilonMatrix_);
        nonConstThis->implicitSolver_->SetRelTol(1e-12);
        nonConstThis->implicitSolver_->SetMaxIter(1000);
        nonConstThis->implicitSolver_->SetPrintLevel(0);
    }
    
    mfem::Vector rhs(dEdt);
    dEdt = 0.0;
    implicitSolver_->Mult(rhs, dEdt);
    
    boundaries_->applyDirichletBC(dEdt);
#endif
}

void PhysicsMaxwellTimeDomain::implicitSolve(double dt, const mfem::Vector& B, 
                                            mfem::Vector& dEdt)
{
    // Implicit solve for lossy case
    // From: ε ∂E/∂t = ∇×(μ⁻¹B) - σE
    // Backward Euler: (M_ε + dt*M_σ) * E^{n+1} = M_ε * E^n + dt * ∇×(μ⁻¹B^{n+1})
    // So: (M_ε + dt*M_σ) * dE/dt = ∇×(μ⁻¹B) - σE^n
    
#ifdef MFEM_USE_MPI
    // Compute RHS: ∇×(μ⁻¹B)
    weakCurlMatrix_->Mult(B, *rhsVec_);
    
    // Subtract source term if present
    if (sources_->hasSource())
    {
        mfem::HypreParVector* jDual = sources_->getDualForm()->ParallelAssemble();
        *rhsVec_ -= *jDual;
        delete jDual;
    }
    
    // Subtract σE^n from RHS
    if (isLossy_ && lossMatrix_ != nullptr)
    {
        mfem::HypreParVector temp(hCurlFESpace_);
        lossMatrix_->Mult(*eVec_, temp);
        *rhsVec_ -= temp;
    }
    
    // Assemble system matrix: A = M_ε + dt*M_σ
    mfem::HypreParMatrix* systemMatrix = nullptr;
    
    if (isLossy_ && lossMatrix_ != nullptr)
    {
        // A = M_ε + dt*M_σ
        systemMatrix = mfem::Add(1.0, *massEpsilonMatrix_, dt, *lossMatrix_);
    }
    else
    {
        // No losses, just use M_ε
        systemMatrix = massEpsilonMatrix_;
    }
    
    // Setup or update solver
    if (implicitSolver_ == nullptr)
    {
        implicitSolver_ = new mfem::CGSolver(pmesh_->GetComm());
        preconditioner_ = new mfem::HypreSmoother();
        preconditioner_->SetType(mfem::HypreSmoother::Jacobi);
        implicitSolver_->SetRelTol(1e-12);
        implicitSolver_->SetMaxIter(1000);
        implicitSolver_->SetPrintLevel(0);
    }
    
    // Update operator and preconditioner
    implicitSolver_->SetOperator(*systemMatrix);
    if (systemMatrix != massEpsilonMatrix_)
    {
        preconditioner_->SetOperator(*systemMatrix);
    }
    else
    {
        preconditioner_->SetOperator(*massEpsilonMatrix_);
    }
    implicitSolver_->SetPreconditioner(*preconditioner_);
    
    // Solve for dE/dt
    dEdt = 0.0;
    implicitSolver_->Mult(*rhsVec_, dEdt);
    
    // Apply boundary conditions
    boundaries_->applyDirichletBC(dEdt);
    
    // Clean up if we created a new matrix
    if (systemMatrix != massEpsilonMatrix_)
    {
        delete systemMatrix;
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
    
    // Subtract σE^n
    if (isLossy_ && lossMatrix_ != nullptr)
    {
        mfem::Vector temp(hCurlFESpace_->GetTrueVSize());
        lossMatrix_->Mult(*eField_, temp);
        dEdt -= temp;
    }
    
    // Assemble system matrix
    mfem::SparseMatrix* systemMatrix = nullptr;
    
    if (isLossy_ && lossMatrix_ != nullptr)
    {
        systemMatrix = mfem::Add(1.0, *massEpsilonMatrix_, dt, *lossMatrix_);
    }
    else
    {
        systemMatrix = massEpsilonMatrix_;
    }
    
    // Setup or update solver
    if (implicitSolver_ == nullptr)
    {
        implicitSolver_ = new mfem::CGSolver();
        implicitSolver_->SetRelTol(1e-12);
        implicitSolver_->SetMaxIter(1000);
        implicitSolver_->SetPrintLevel(0);
    }
    
    implicitSolver_->SetOperator(*systemMatrix);
    
    // Solve
    mfem::Vector rhs(dEdt);
    dEdt = 0.0;
    implicitSolver_->Mult(rhs, dEdt);
    
    boundaries_->applyDirichletBC(dEdt);
    
    // Clean up
    if (systemMatrix != massEpsilonMatrix_)
    {
        delete systemMatrix;
    }
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
    // Compute electromagnetic energy: E = 0.5 * (E^T M_ε E + B^T M_μ^{-1} B)
    
    double energyE = 0.0;
    double energyB = 0.0;
    
#ifdef MFEM_USE_MPI
    // Electric field energy: 0.5 * E^T M_ε E
    mfem::HypreParVector tempE(hCurlFESpace_);
    massEpsilonMatrix_->Mult(*eVec_, tempE);
    energyE = 0.5 * mfem::InnerProduct(*eVec_, tempE);
    
    // Magnetic field energy: 0.5 * B^T M_μ^{-1} B
    mfem::HypreParVector tempB(hDivFESpace_);
    massMuInvMatrix_->Mult(*bVec_, tempB);
    energyB = 0.5 * mfem::InnerProduct(*bVec_, tempB);
    
    // InnerProduct already returns global sum, no need for additional Allreduce
    return energyE + energyB;
#else
    // Serial version
    mfem::Vector tempE(hCurlFESpace_->GetTrueVSize());
    massEpsilonMatrix_->Mult(*eField_, tempE);
    energyE = 0.5 * (*eField_) * tempE;
    
    mfem::Vector tempB(hDivFESpace_->GetTrueVSize());
    massMuInvMatrix_->Mult(*bField_, tempB);
    energyB = 0.5 * (*bField_) * tempB;
    
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

void PhysicsMaxwellTimeDomain::updateSourceTime(double time)
{
    if (sources_->hasSource())
    {
        sources_->updateTime(time);
    }
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
