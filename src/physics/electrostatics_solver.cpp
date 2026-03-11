#include "electrostatics_solver.hpp"
#include "logger.hpp"

#include <cmath>

namespace mpfem {

namespace {
const double MIN_RESISTIVITY = 1e-16;
const double DEFAULT_REFERENCE_TEMPERATURE = 298;
} // namespace

// --- ConductivityCoefficient ---

void ConductivityCoefficient::setMaterialFields(const mfem::Vector& rho0,
                                                 const mfem::Vector& alpha,
                                                 const mfem::Vector& tref,
                                                 const mfem::Vector& sigma0)
{
    rho0_ = rho0;
    alpha_ = alpha;
    tref_ = tref;
    sigma0_ = sigma0;
}

void ConductivityCoefficient::setTemperatureField(const mfem::GridFunction* temperature)
{
    temperature_ = temperature;
}

double ConductivityCoefficient::Eval(mfem::ElementTransformation& tr,
                                     const mfem::IntegrationPoint& ip)
{
    const int attr = tr.Attribute;
    if (attr <= 0 || attr > sigma0_.Size()) {
        Logger::log(LogLevel::Error, "ConductivityCoefficient: invalid attribute " 
            + std::to_string(attr) + " (max=" + std::to_string(sigma0_.Size()) + ")");
        return 1.0;
    }

    const double rho0 = rho0_(attr - 1);
    if (rho0 > 0.0) {
        // Linear resistivity model: rho = rho0 * (1 + alpha * (T - Tref))
        const double alpha = alpha_(attr - 1);
        const double tref = tref_(attr - 1);
        double temp = tref;
        if (temperature_) {
            temp = temperature_->GetValue(tr, ip);
        }
        const double rho = rho0 * (1.0 + alpha * (temp - tref));
        Check(std::isfinite(rho) && rho > 0.0,
              "ConductivityCoefficient: resistivity is invalid (attr=" +
              std::to_string(attr) + ", rho=" + std::to_string(rho) + ")");
        return 1.0 / rho;
    }

    const double sigma = sigma0_(attr - 1);
    Check(std::isfinite(sigma) && sigma > 0.0,
          "ConductivityCoefficient: conductivity is invalid (attr=" +
          std::to_string(attr) + ", value=" + std::to_string(sigma) + ")");
    return sigma;
}

// --- ElectrostaticsSolver ---

void ElectrostaticsSolver::setOrder(int order)
{
    order_ = order > 0 ? order : 1;
}

void ElectrostaticsSolver::setSolver(std::unique_ptr<LinearSolverStrategy> solver)
{
    solver_ = std::move(solver);
}

void ElectrostaticsSolver::fillMaterialVectors(const PhysicsMaterialDatabase& materials)
{
    const int maxAttr = mesh_->attributes.Max();
    rho0_.SetSize(maxAttr);
    alpha_.SetSize(maxAttr);
    tref_.SetSize(maxAttr);
    sigma0_.SetSize(maxAttr);

    rho0_ = 0.0;
    alpha_ = 0.0;
    tref_ = DEFAULT_REFERENCE_TEMPERATURE;
    sigma0_ = 1.0;

    for (int attr = 1; attr <= maxAttr; ++attr) {
        auto tagIt = problemModel_->domainMaterialTag.find(attr);
        if (tagIt == problemModel_->domainMaterialTag.end()) continue;

        auto propIt = materials.byTag.find(tagIt->second);
        if (propIt == materials.byTag.end()) continue;

        const auto& prop = propIt->second;
        rho0_(attr - 1) = prop.rho0;
        alpha_(attr - 1) = prop.alpha;
        tref_(attr - 1) = prop.tref;
        sigma0_(attr - 1) = prop.electricConductivity > 0.0 ? prop.electricConductivity : 1.0;
    }
}

void ElectrostaticsSolver::buildBoundaryMarkers()
{
    const int maxBdrAttr = mesh_->bdr_attributes.Max();
    essentialBdr_.SetSize(maxBdrAttr);
    essentialBdr_ = 0;
    dirichletValues_.SetSize(maxBdrAttr);
    dirichletValues_ = 0.0;

    auto fieldIt = problemModel_->boundaries.find(FieldKind::ElectricPotential);
    if (fieldIt == problemModel_->boundaries.end()) return;

    for (const auto& [id, params] : fieldIt->second) {
        if (params.kind == "voltage" && id > 0 && id <= maxBdrAttr) {
            essentialBdr_[id - 1] = 1;
            auto valIt = params.values.find("value");
            dirichletValues_(id - 1) = (valIt != params.values.end()) ? valIt->second : 0.0;
        }
    }
}

void ElectrostaticsSolver::initialize(mfem::Mesh& mesh,
                                      const PhysicsProblemModel& problemModel,
                                      const PhysicsMaterialDatabase& materials)
{
    mesh_ = &mesh;
    problemModel_ = &problemModel;

    // Create finite element space
    fec_ = std::make_unique<mfem::H1_FECollection>(order_, mesh.Dimension());
    space_ = std::make_unique<mfem::FiniteElementSpace>(&mesh, fec_.get());

    // Initialize solution
    potential_ = std::make_unique<mfem::GridFunction>(space_.get());
    *potential_ = 0.0;

    // Material properties
    fillMaterialVectors(materials);
    conductivity_.setMaterialFields(rho0_, alpha_, tref_, sigma0_);

    // Boundary markers
    buildBoundaryMarkers();

    // Cache essential true DOFs (doesn't change between iterations)
    space_->GetEssentialTrueDofs(essentialBdr_, essTdof_);
}

void ElectrostaticsSolver::applyBoundaryConditions()
{
    mfem::PWConstCoefficient boundaryCoef(dirichletValues_);
    potential_->ProjectBdrCoefficient(boundaryCoef, essentialBdr_);
}

void ElectrostaticsSolver::assemble()
{
    // Create fresh forms each iteration (MFEM accumulates values on re-assembly)
    // The forms themselves are cheap to create; the expensive part is element iteration
    // which must happen anyway for coefficient evaluation.
    aForm_ = std::make_unique<mfem::BilinearForm>(space_.get());
    aForm_->AddDomainIntegrator(new mfem::DiffusionIntegrator(conductivity_));
    aForm_->Assemble();

    bForm_ = std::make_unique<mfem::LinearForm>(space_.get());
    bForm_->Assemble();
}

void ElectrostaticsSolver::solve()
{
    Check(solver_, "No linear solver set for ElectrostaticsSolver");

    mfem::SparseMatrix mat;
    mfem::Vector x, b;
    aForm_->FormLinearSystem(essTdof_, *potential_, *bForm_, mat, x, b);
    x = 0.0;

    solver_->solve(mat, x, b);

    aForm_->RecoverFEMSolution(x, *bForm_, *potential_);
}

mfem::GridFunction& ElectrostaticsSolver::getField() { return *potential_; }
const mfem::GridFunction& ElectrostaticsSolver::getField() const { return *potential_; }
FieldKind ElectrostaticsSolver::getFieldKind() const { return FieldKind::ElectricPotential; }
mfem::FiniteElementSpace& ElectrostaticsSolver::getSpace() { return *space_; }
const mfem::FiniteElementSpace& ElectrostaticsSolver::getSpace() const { return *space_; }

void ElectrostaticsSolver::setTemperatureField(const mfem::GridFunction* temperature)
{
    conductivity_.setTemperatureField(temperature);
}

mfem::Coefficient* ElectrostaticsSolver::getConductivityCoefficient()
{
    return &conductivity_;
}

} // namespace mpfem