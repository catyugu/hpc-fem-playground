/**
 * @file physics/physics_joule_heating.cpp
 * @brief Implementation of coupled Joule heating physics with nonlinear coupling
 */

#include "hpcfem/physics/physics_joule_heating.hpp"
#include "joule_heating_coefficient.hpp"

namespace hpcfem
{

#ifdef MFEM_USE_MPI
JouleHeatingPhysics::JouleHeatingPhysics(mfem::ParMesh* pmesh,
                                         int polynomialOrder,
                                         double electricalConductivity,
                                         double thermalConductivity)
    : sigma_(electricalConductivity),
      kappa_(thermalConductivity),
      matrixK_e_(nullptr),
      matrixK_t_(nullptr),
      matrixC_(nullptr)
{
    int dim = pmesh->Dimension();
    fecElectric_ = new mfem::H1_FECollection(polynomialOrder, dim);
    fecThermal_ = new mfem::H1_FECollection(polynomialOrder, dim);
    fespaceElectric_ = new mfem::ParFiniteElementSpace(pmesh, fecElectric_);
    fespaceThermal_ = new mfem::ParFiniteElementSpace(pmesh, fecThermal_);
    blockOffsets_.SetSize(3);
    blockOffsets_[0] = 0;
    blockOffsets_[1] = fespaceElectric_->GetTrueVSize();
    blockOffsets_[2] = blockOffsets_[1] + fespaceThermal_->GetTrueVSize();
    sigmaCoeff_ = new mfem::ConstantCoefficient(sigma_);
    kappaCoeff_ = new mfem::ConstantCoefficient(kappa_);
    bilinearElectric_ = new mfem::ParBilinearForm(fespaceElectric_);
    bilinearElectric_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*sigmaCoeff_));
    bilinearThermal_ = new mfem::ParBilinearForm(fespaceThermal_);
    bilinearThermal_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*kappaCoeff_));
    linearElectric_ = new mfem::ParLinearForm(fespaceElectric_);
    linearThermal_ = new mfem::ParLinearForm(fespaceThermal_);
}
#else
JouleHeatingPhysics::JouleHeatingPhysics(mfem::Mesh* mesh,
                                         int polynomialOrder,
                                         double electricalConductivity,
                                         double thermalConductivity)
    : sigma_(electricalConductivity),
      kappa_(thermalConductivity),
      matrixK_e_(nullptr),
      matrixK_t_(nullptr),
      matrixC_(nullptr)
{
    int dim = mesh->Dimension();
    fecElectric_ = new mfem::H1_FECollection(polynomialOrder, dim);
    fecThermal_ = new mfem::H1_FECollection(polynomialOrder, dim);
    fespaceElectric_ = new mfem::FiniteElementSpace(mesh, fecElectric_);
    fespaceThermal_ = new mfem::FiniteElementSpace(mesh, fecThermal_);
    blockOffsets_.SetSize(3);
    blockOffsets_[0] = 0;
    blockOffsets_[1] = fespaceElectric_->GetTrueVSize();
    blockOffsets_[2] = blockOffsets_[1] + fespaceThermal_->GetTrueVSize();
    sigmaCoeff_ = new mfem::ConstantCoefficient(sigma_);
    kappaCoeff_ = new mfem::ConstantCoefficient(kappa_);
    bilinearElectric_ = new mfem::BilinearForm(fespaceElectric_);
    bilinearElectric_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*sigmaCoeff_));
    bilinearThermal_ = new mfem::BilinearForm(fespaceThermal_);
    bilinearThermal_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*kappaCoeff_));
    linearElectric_ = new mfem::LinearForm(fespaceElectric_);
    linearThermal_ = new mfem::LinearForm(fespaceThermal_);
}
#endif

JouleHeatingPhysics::~JouleHeatingPhysics()
{
    delete matrixC_;
    delete matrixK_t_;
    delete matrixK_e_;
    delete linearThermal_;
    delete linearElectric_;
    delete bilinearThermal_;
    delete bilinearElectric_;
    delete kappaCoeff_;
    delete sigmaCoeff_;
    delete fespaceThermal_;
    delete fespaceElectric_;
    delete fecThermal_;
    delete fecElectric_;
}

#ifdef MFEM_USE_MPI
void JouleHeatingPhysics::assemble(mfem::BlockOperator& blockOperator,
                                   mfem::BlockVector& blockRHS,
                                   mfem::BlockVector& blockSolution,
                                   mfem::Array<int>& essTdofElectric,
                                   mfem::Array<int>& essTdofThermal,
                                   mfem::ParGridFunction* voltageGF)
{
    delete matrixK_e_;
    delete matrixK_t_;
    matrixK_e_ = nullptr;
    matrixK_t_ = nullptr;
    delete bilinearElectric_;
    delete bilinearThermal_;
    delete linearElectric_;
    delete linearThermal_;
    bilinearElectric_ = new mfem::ParBilinearForm(fespaceElectric_);
    bilinearElectric_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*sigmaCoeff_));
    bilinearThermal_ = new mfem::ParBilinearForm(fespaceThermal_);
    bilinearThermal_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*kappaCoeff_));
    bilinearElectric_->Assemble();
    bilinearElectric_->Finalize();
    bilinearThermal_->Assemble();
    bilinearThermal_->Finalize();
    linearElectric_ = new mfem::ParLinearForm(fespaceElectric_);
    linearElectric_->Assemble();
    linearThermal_ = new mfem::ParLinearForm(fespaceThermal_);
    if (voltageGF != nullptr) {
        JouleHeatingCoefficient jouleCoeff(voltageGF, sigma_);
        linearThermal_->AddDomainIntegrator(
            new mfem::DomainLFIntegrator(jouleCoeff));
    }
    linearThermal_->Assemble();
    mfem::ParGridFunction electricBC(fespaceElectric_);
    mfem::ParGridFunction thermalBC(fespaceThermal_);
    electricBC = 0.0;
    thermalBC = 0.0;
    mfem::Array<int> essBdrElectric(fespaceElectric_->GetParMesh()->bdr_attributes.Max());
    essBdrElectric = 1;
    fespaceElectric_->GetEssentialTrueDofs(essBdrElectric, essTdofElectric);
    mfem::Array<int> essBdrThermal(fespaceThermal_->GetParMesh()->bdr_attributes.Max());
    essBdrThermal = 1;
    fespaceThermal_->GetEssentialTrueDofs(essBdrThermal, essTdofThermal);
    mfem::Vector rhsElectric, rhsThermal;
    mfem::Vector solElectric, solThermal;
    mfem::HypreParMatrix A_electric, A_thermal;
    bilinearElectric_->FormLinearSystem(essTdofElectric, electricBC, *linearElectric_,
                                       A_electric, solElectric, rhsElectric);
    bilinearThermal_->FormLinearSystem(essTdofThermal, thermalBC, *linearThermal_,
                                      A_thermal, solThermal, rhsThermal);
    matrixK_e_ = new mfem::HypreParMatrix(A_electric);
    matrixK_t_ = new mfem::HypreParMatrix(A_thermal);
    blockOperator.SetBlock(0, 0, matrixK_e_);
    blockOperator.SetBlock(1, 1, matrixK_t_);
    blockRHS.GetBlock(0) = rhsElectric;
    blockRHS.GetBlock(1) = rhsThermal;
    blockSolution.GetBlock(0) = solElectric;
    blockSolution.GetBlock(1) = solThermal;
}
#else
void JouleHeatingPhysics::assemble(mfem::BlockOperator& blockOperator,
                                   mfem::BlockVector& blockRHS,
                                   mfem::BlockVector& blockSolution,
                                   mfem::Array<int>& essTdofElectric,
                                   mfem::Array<int>& essTdofThermal,
                                   mfem::GridFunction* voltageGF)
{
    delete matrixK_e_;
    delete matrixK_t_;
    matrixK_e_ = nullptr;
    matrixK_t_ = nullptr;
    delete bilinearElectric_;
    delete bilinearThermal_;
    delete linearElectric_;
    delete linearThermal_;
    bilinearElectric_ = new mfem::BilinearForm(fespaceElectric_);
    bilinearElectric_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*sigmaCoeff_));
    bilinearThermal_ = new mfem::BilinearForm(fespaceThermal_);
    bilinearThermal_->AddDomainIntegrator(new mfem::DiffusionIntegrator(*kappaCoeff_));
    bilinearElectric_->Assemble();
    bilinearElectric_->Finalize();
    bilinearThermal_->Assemble();
    bilinearThermal_->Finalize();
    linearElectric_ = new mfem::LinearForm(fespaceElectric_);
    linearElectric_->Assemble();
    linearThermal_ = new mfem::LinearForm(fespaceThermal_);
    if (voltageGF != nullptr) {
        JouleHeatingCoefficient jouleCoeff(voltageGF, sigma_);
        linearThermal_->AddDomainIntegrator(
            new mfem::DomainLFIntegrator(jouleCoeff));
    }
    linearThermal_->Assemble();
    mfem::GridFunction electricBC(fespaceElectric_);
    mfem::GridFunction thermalBC(fespaceThermal_);
    electricBC = 0.0;
    thermalBC = 0.0;
    mfem::Array<int> essBdrElectric(fespaceElectric_->GetMesh()->bdr_attributes.Max());
    essBdrElectric = 1;
    fespaceElectric_->GetEssentialTrueDofs(essBdrElectric, essTdofElectric);
    mfem::Array<int> essBdrThermal(fespaceThermal_->GetMesh()->bdr_attributes.Max());
    essBdrThermal = 1;
    fespaceThermal_->GetEssentialTrueDofs(essBdrThermal, essTdofThermal);
    mfem::Vector rhsElectric, rhsThermal;
    mfem::Vector solElectric, solThermal;
    mfem::SparseMatrix A_electric, A_thermal;
    bilinearElectric_->FormLinearSystem(essTdofElectric, electricBC, *linearElectric_,
                                       A_electric, solElectric, rhsElectric);
    bilinearThermal_->FormLinearSystem(essTdofThermal, thermalBC, *linearThermal_,
                                      A_thermal, solThermal, rhsThermal);
    matrixK_e_ = new mfem::SparseMatrix(A_electric);
    matrixK_t_ = new mfem::SparseMatrix(A_thermal);
    blockOperator.SetBlock(0, 0, matrixK_e_);
    blockOperator.SetBlock(1, 1, matrixK_t_);
    blockRHS.GetBlock(0) = rhsElectric;
    blockRHS.GetBlock(1) = rhsThermal;
    blockSolution.GetBlock(0) = solElectric;
    blockSolution.GetBlock(1) = solThermal;
}
#endif

#ifdef MFEM_USE_MPI
mfem::ParFiniteElementSpace* JouleHeatingPhysics::getElectricSpace()
{
    return fespaceElectric_;
}

mfem::ParFiniteElementSpace* JouleHeatingPhysics::getThermalSpace()
{
    return fespaceThermal_;
}
#else
mfem::FiniteElementSpace* JouleHeatingPhysics::getElectricSpace()
{
    return fespaceElectric_;
}

mfem::FiniteElementSpace* JouleHeatingPhysics::getThermalSpace()
{
    return fespaceThermal_;
}
#endif

const mfem::Array<int>& JouleHeatingPhysics::getBlockOffsets() const
{
    return blockOffsets_;
}

} // namespace hpcfem
