#include "solid_mechanics_solver.hpp"
#include "logger.hpp"

#include <cmath>

namespace mpfem {

namespace {
const double DEFAULT_REFERENCE_TEMPERATURE = 293.15;
} // namespace

// --- ThermalExpansionLoadCoefficient ---

void ThermalExpansionLoadCoefficient::Eval(mfem::Vector& V,
                                           mfem::ElementTransformation& tr,
                                           const mfem::IntegrationPoint& ip)
{
    V.SetSize(dimension_ * dimension_);
    V = 0.0;

    if (!temperature_ || !alpha_ || !lambda_ || !mu_) return;

    const double temp = temperature_->GetValue(tr, ip);
    const double alpha = alpha_->Eval(tr, ip);
    const double lambda = lambda_->Eval(tr, ip);
    const double mu = mu_->Eval(tr, ip);
    const double deltaT = temp - referenceTemperature_;

    const double thermalModulus = (3.0 * lambda + 2.0 * mu) * alpha * deltaT;

    for (int i = 0; i < dimension_; ++i) {
        V(i * dimension_ + i) = thermalModulus;
    }
}

// --- SolidMechanicsSolver ---

void SolidMechanicsSolver::setOrder(int order)
{
    order_ = order > 0 ? order : 1;
}

void SolidMechanicsSolver::setSolver(std::unique_ptr<LinearSolverStrategy> solver)
{
    solver_ = std::move(solver);
}

void SolidMechanicsSolver::fillMaterialVectors(const PhysicsMaterialDatabase& materials)
{
    const int maxAttr = mesh_->attributes.Max();
    young_.SetSize(maxAttr);
    poisson_.SetSize(maxAttr);
    thermalExpansion_.SetSize(maxAttr);

    young_ = 1.0;
    poisson_ = 0.3;
    thermalExpansion_ = 1.0e-5;

    for (int attr = 1; attr <= maxAttr; ++attr) {
        auto tagIt = problemModel_->domainMaterialTag.find(attr);
        if (tagIt == problemModel_->domainMaterialTag.end()) continue;

        auto propIt = materials.byTag.find(tagIt->second);
        if (propIt == materials.byTag.end()) continue;

        const auto& prop = propIt->second;
        young_(attr - 1) = prop.youngModulus > 0.0 ? prop.youngModulus : 1.0;
        poisson_(attr - 1) = prop.poissonRatio;
        thermalExpansion_(attr - 1) = prop.thermalExpansion > 0.0 ? prop.thermalExpansion : 1.0e-5;
    }
}

void SolidMechanicsSolver::computeLameParameters()
{
    const int n = young_.Size();
    lambdaValues_.SetSize(n);
    muValues_.SetSize(n);

    for (int i = 0; i < n; ++i) {
        const double e = young_(i);
        const double nu = poisson_(i);
        muValues_(i) = e / (2.0 * (1.0 + nu));
        lambdaValues_(i) = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    }

    lambda_.setValues(lambdaValues_);
    mu_.setValues(muValues_);
    alpha_.setValues(thermalExpansion_);
}

void SolidMechanicsSolver::buildBoundaryMarkers()
{
    const int maxBdrAttr = mesh_->bdr_attributes.Max();
    fixedBdr_.SetSize(maxBdrAttr);
    fixedBdr_ = 0;

    auto fieldIt = problemModel_->boundaries.find(FieldKind::Displacement);
    if (fieldIt == problemModel_->boundaries.end()) return;

    for (const auto& [id, params] : fieldIt->second) {
        if (params.kind == "fixed_constraint" && id > 0 && id <= maxBdrAttr) {
            fixedBdr_[id - 1] = 1;
        }
    }
}

void SolidMechanicsSolver::parseCouplingTerms()
{
    hasThermalExpansion_ = false;
    for (const auto& coupling : problemModel_->couplings) {
        if (coupling == CouplingKind::ThermalExpansion) {
            hasThermalExpansion_ = true;
            break;
        }
    }
}

void SolidMechanicsSolver::initialize(mfem::Mesh& mesh,
                                      const PhysicsProblemModel& problemModel,
                                      const PhysicsMaterialDatabase& materials)
{
    mesh_ = &mesh;
    problemModel_ = &problemModel;

    // Create vector finite element space
    fec_ = std::make_unique<mfem::H1_FECollection>(order_, mesh.Dimension());
    space_ = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec_.get(), mesh.Dimension());

    // Initialize solution
    displacement_ = std::make_unique<mfem::GridFunction>(space_.get());
    *displacement_ = 0.0;

    // Material properties
    fillMaterialVectors(materials);
    computeLameParameters();

    // Boundary conditions
    buildBoundaryMarkers();

    // Coupling terms
    parseCouplingTerms();

    // Cache essential true DOFs
    space_->GetEssentialTrueDofs(fixedBdr_, essTdof_);

    // Zero displacement for fixed boundaries
    zeroDisp_.SetSize(mesh.Dimension());
    zeroDisp_ = 0.0;
    zeroCoef_ = std::make_unique<mfem::VectorConstantCoefficient>(zeroDisp_);

    // Thermal expansion load coefficient
    if (hasThermalExpansion_) {
        thermalExpansionLoad_.emplace(mesh.Dimension());
        thermalExpansionLoad_->setAlphaCoefficient(&alpha_);
        thermalExpansionLoad_->setLambdaCoefficient(&lambda_);
        thermalExpansionLoad_->setMuCoefficient(&mu_);
        thermalExpansionLoad_->setReferenceTemperature(referenceTemperature_);
    }
}

void SolidMechanicsSolver::applyBoundaryConditions()
{
    displacement_->ProjectBdrCoefficient(*zeroCoef_, fixedBdr_);
}

void SolidMechanicsSolver::assemble()
{
    // Create fresh forms each iteration
    aForm_ = std::make_unique<mfem::BilinearForm>(space_.get());
    aForm_->AddDomainIntegrator(new mfem::ElasticityIntegrator(lambda_, mu_));
    aForm_->Assemble();

    bForm_ = std::make_unique<mfem::LinearForm>(space_.get());
    if (hasThermalExpansion_ && thermalExpansionLoad_) {
        bForm_->AddDomainIntegrator(new mfem::VectorDomainLFGradIntegrator(*thermalExpansionLoad_));
    }
    bForm_->Assemble();
}

void SolidMechanicsSolver::solve()
{
    Check(solver_, "No linear solver set for SolidMechanicsSolver");

    mfem::SparseMatrix mat;
    mfem::Vector x, b;
    aForm_->FormLinearSystem(essTdof_, *displacement_, *bForm_, mat, x, b);
    x = 0.0;

    solver_->solve(mat, x, b);

    aForm_->RecoverFEMSolution(x, *bForm_, *displacement_);
}

mfem::GridFunction& SolidMechanicsSolver::getField() { return *displacement_; }
const mfem::GridFunction& SolidMechanicsSolver::getField() const { return *displacement_; }
FieldKind SolidMechanicsSolver::getFieldKind() const { return FieldKind::Displacement; }
mfem::FiniteElementSpace& SolidMechanicsSolver::getSpace() { return *space_; }
const mfem::FiniteElementSpace& SolidMechanicsSolver::getSpace() const { return *space_; }

void SolidMechanicsSolver::setTemperatureField(const mfem::GridFunction* temperature)
{
    if (thermalExpansionLoad_) {
        thermalExpansionLoad_->setTemperatureField(temperature);
    }
}

mfem::Coefficient* SolidMechanicsSolver::getLambdaCoefficient() { return &lambda_; }
mfem::Coefficient* SolidMechanicsSolver::getMuCoefficient() { return &mu_; }
mfem::Coefficient* SolidMechanicsSolver::getThermalExpansionCoefficient() { return &alpha_; }

} // namespace mpfem