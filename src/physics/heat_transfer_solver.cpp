#include "heat_transfer_solver.hpp"
#include "logger.hpp"

#include <cmath>

namespace mpfem {

namespace {
const double DEFAULT_REFERENCE_TEMPERATURE = 293.15;
} // namespace

// --- AttributeCoefficient ---

double AttributeCoefficient::Eval(mfem::ElementTransformation& tr,
                                  const mfem::IntegrationPoint& ip)
{
    (void)ip;
    const int attr = tr.Attribute;
    if (attr <= 0 || attr > values_.Size()) return 1.0;
    return values_(attr - 1);
}

// --- JouleHeatCoefficient ---

double JouleHeatCoefficient::Eval(mfem::ElementTransformation& tr,
                                  const mfem::IntegrationPoint& ip)
{
    Check(potential_ != nullptr, "JouleHeatCoefficient: potential field is null");
    Check(conductivity_ != nullptr, "JouleHeatCoefficient: conductivity coefficient is null");

    mfem::Vector grad;
    potential_->GetGradient(tr, grad);

    for (int i = 0; i < grad.Size(); ++i) {
        Check(std::isfinite(grad(i)),
              "JouleHeatCoefficient: gradient component " + std::to_string(i) +
              " is not finite (value: " + std::to_string(grad(i)) + ")");
    }

    const double sigma = conductivity_->Eval(tr, ip);
    Check(std::isfinite(sigma),
          "JouleHeatCoefficient: conductivity is not finite (value: " + std::to_string(sigma) + ")");

    const double norm2 = grad * grad;
    const double result = sigma * norm2;
    Check(std::isfinite(result) && result >= 0.0,
          "JouleHeatCoefficient: result is invalid (sigma=" + std::to_string(sigma) +
          ", norm2=" + std::to_string(norm2) + ", result=" + std::to_string(result) + ")");

    return result;
}

// --- HeatTransferSolver ---

void HeatTransferSolver::setOrder(int order)
{
    order_ = order > 0 ? order : 1;
}

void HeatTransferSolver::setSolver(std::unique_ptr<LinearSolverStrategy> solver)
{
    solver_ = std::move(solver);
}

void HeatTransferSolver::fillMaterialVectors(const PhysicsMaterialDatabase& materials)
{
    const int maxAttr = mesh_->attributes.Max();
    thermalK_.SetSize(maxAttr);
    thermalK_ = 1.0;

    for (int attr = 1; attr <= maxAttr; ++attr) {
        auto tagIt = problemModel_->domainMaterialTag.find(attr);
        if (tagIt == problemModel_->domainMaterialTag.end()) continue;

        auto propIt = materials.byTag.find(tagIt->second);
        if (propIt == materials.byTag.end()) continue;

        const auto& prop = propIt->second;
        thermalK_(attr - 1) = prop.thermalConductivity > 0.0 ? prop.thermalConductivity : 1.0;
    }
}

void HeatTransferSolver::buildConvectionBoundaryMarkers()
{
    const int maxBdrAttr = mesh_->bdr_attributes.Max();
    convectionBdr_.SetSize(maxBdrAttr);
    convectionBdr_ = 0;
    convectionH_.SetSize(maxBdrAttr);
    convectionH_ = 0.0;
    convectionTinf_.SetSize(maxBdrAttr);
    convectionTinf_ = DEFAULT_REFERENCE_TEMPERATURE;

    auto fieldIt = problemModel_->boundaries.find(FieldKind::Temperature);
    if (fieldIt == problemModel_->boundaries.end()) return;

    for (const auto& [id, params] : fieldIt->second) {
        if (params.kind == "convection" && id > 0 && id <= maxBdrAttr) {
            convectionBdr_[id - 1] = 1;
            auto hIt = params.values.find("h");
            auto tinfIt = params.values.find("T_inf");
            convectionH_(id - 1) = (hIt != params.values.end()) ? hIt->second : 0.0;
            convectionTinf_(id - 1) = (tinfIt != params.values.end()) 
                ? tinfIt->second : DEFAULT_REFERENCE_TEMPERATURE;
        }
    }
}

void HeatTransferSolver::parseSourceTerms()
{
    hasJouleHeating_ = false;
    for (const auto& source : problemModel_->sources) {
        if (source.field == FieldKind::Temperature &&
            source.coupled && source.couplingKind == CouplingKind::JouleHeating) {
            hasJouleHeating_ = true;
            break;
        }
    }
}

void HeatTransferSolver::initialize(mfem::Mesh& mesh,
                                    const PhysicsProblemModel& problemModel,
                                    const PhysicsMaterialDatabase& materials)
{
    mesh_ = &mesh;
    problemModel_ = &problemModel;

    // Create finite element space
    fec_ = std::make_unique<mfem::H1_FECollection>(order_, mesh.Dimension());
    space_ = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec_.get());

    // Initialize solution
    temperature_ = std::make_unique<mfem::GridFunction>(space_.get());
    *temperature_ = DEFAULT_REFERENCE_TEMPERATURE;

    // Material properties
    fillMaterialVectors(materials);
    thermalConductivity_.setValues(thermalK_);

    // Boundary conditions
    buildConvectionBoundaryMarkers();

    // Source terms
    parseSourceTerms();

    // Create convection coefficients (constant)
    if (convectionBdr_.Max() > 0) {
        hCoef_ = std::make_unique<mfem::PWConstCoefficient>(convectionH_);
        tinfCoef_ = std::make_unique<mfem::PWConstCoefficient>(convectionTinf_);
        htinfCoef_ = std::make_unique<mfem::ProductCoefficient>(*hCoef_, *tinfCoef_);
    }

    // No essential DOFs for heat transfer (only Neumann/Robin BCs)
    essTdof_.SetSize(0);
}

void HeatTransferSolver::applyBoundaryConditions()
{
    // Heat transfer has no essential BCs (only convection)
}

void HeatTransferSolver::assemble()
{
    // Create fresh forms each iteration (MFEM accumulates values on re-assembly)
    aForm_ = std::make_unique<mfem::BilinearForm>(space_.get());
    aForm_->AddDomainIntegrator(new mfem::DiffusionIntegrator(thermalConductivity_));
    if (convectionBdr_.Max() > 0) {
        aForm_->AddBoundaryIntegrator(new mfem::MassIntegrator(*hCoef_), convectionBdr_);
    }
    aForm_->Assemble();

    bForm_ = std::make_unique<mfem::LinearForm>(space_.get());
    if (hasJouleHeating_) {
        bForm_->AddDomainIntegrator(new mfem::DomainLFIntegrator(jouleHeatSource_));
    }
    if (convectionBdr_.Max() > 0) {
        bForm_->AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(*htinfCoef_), convectionBdr_);
    }
    bForm_->Assemble();
}

void HeatTransferSolver::solve()
{
    Check(solver_, "No linear solver set for HeatTransferSolver");

    mfem::SparseMatrix mat;
    mfem::Vector x, b;
    aForm_->FormLinearSystem(essTdof_, *temperature_, *bForm_, mat, x, b);
    x = 0.0;

    solver_->solve(mat, x, b);

    aForm_->RecoverFEMSolution(x, *bForm_, *temperature_);
}

mfem::GridFunction& HeatTransferSolver::getField() { return *temperature_; }
const mfem::GridFunction& HeatTransferSolver::getField() const { return *temperature_; }
FieldKind HeatTransferSolver::getFieldKind() const { return FieldKind::Temperature; }
mfem::FiniteElementSpace& HeatTransferSolver::getSpace() { return *space_; }
const mfem::FiniteElementSpace& HeatTransferSolver::getSpace() const { return *space_; }

void HeatTransferSolver::setPotentialField(const mfem::GridFunction* potential)
{
    jouleHeatSource_.setPotential(potential);
}

void HeatTransferSolver::setConductivityCoefficient(mfem::Coefficient* conductivity)
{
    jouleHeatSource_.setConductivity(conductivity);
}

mfem::Coefficient* HeatTransferSolver::getThermalConductivityCoefficient()
{
    return &thermalConductivity_;
}

} // namespace mpfem