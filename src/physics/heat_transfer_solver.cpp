#include "heat_transfer_solver.hpp"
#include "linear_solver_strategy.hpp"
#include "mfem.hpp"
#include "logger.hpp"

#include <cmath>

namespace mpfem {

namespace {

const double DEFAULT_REFERENCE_TEMPERATURE = 293.15;

/**
 * @brief Attribute-based scalar coefficient.
 */
class AttributeScalarCoefficient : public mfem::Coefficient {
public:
    explicit AttributeScalarCoefficient() : values_() {}

    void setValues(const mfem::Vector& values)
    {
        values_ = values;
    }

    double Eval(mfem::ElementTransformation& transformation,
                const mfem::IntegrationPoint& integrationPoint) override
    {
        (void)integrationPoint;
        const int attribute = transformation.Attribute;
        if (attribute <= 0 || attribute > values_.Size()) {
            return 1.0;
        }
        return values_(attribute - 1);
    }

private:
    mfem::Vector values_;
};

/**
 * @brief Joule heating coefficient: Q = sigma * |grad(V)|^2
 */
class JouleHeatCoefficient : public mfem::Coefficient {
public:
    JouleHeatCoefficient()
        : potential_(nullptr), conductivity_(nullptr) {}

    void setPotential(const mfem::GridFunction* potential)
    {
        potential_ = potential;
    }

    void setConductivity(mfem::Coefficient* conductivity)
    {
        conductivity_ = conductivity;
    }

    double Eval(mfem::ElementTransformation& transformation,
                const mfem::IntegrationPoint& integrationPoint) override
    {
        Check(potential_ != nullptr, "JouleHeatCoefficient: potential field is null");
        Check(conductivity_ != nullptr, "JouleHeatCoefficient: conductivity coefficient is null");

        transformation.SetIntPoint(&integrationPoint);
        mfem::Vector gradient;
        potential_->GetGradient(transformation, gradient);

        // Validate gradient components
        for (int i = 0; i < gradient.Size(); ++i) {
            Check(std::isfinite(gradient(i)),
                  "JouleHeatCoefficient: gradient component " + std::to_string(i) +
                  " is not finite (value: " + std::to_string(gradient(i)) + ")");
        }

        const double sigma = conductivity_->Eval(transformation, integrationPoint);
        Check(std::isfinite(sigma),
              "JouleHeatCoefficient: conductivity is not finite (value: " +
              std::to_string(sigma) + ")");

        const double norm2 = gradient * gradient;
        Check(std::isfinite(norm2),
              "JouleHeatCoefficient: gradient norm squared is not finite (value: " +
              std::to_string(norm2) + ")");

        return sigma * norm2;
    }

private:
    const mfem::GridFunction* potential_;
    mfem::Coefficient* conductivity_;
};

} // namespace

class HeatTransferSolver::Impl {
public:
    int order_ = 1;
    mfem::Mesh* mesh_ = nullptr;
    std::unique_ptr<mfem::H1_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> space_;
    std::unique_ptr<mfem::GridFunction> temperature_;
    std::unique_ptr<LinearSolverStrategy> solver_;

    std::unique_ptr<AttributeScalarCoefficient> thermalConductivity_;
    mfem::Vector thermalK_;

    // Convection boundary conditions
    mfem::Array<int> convectionBdr_;
    mfem::Vector convectionH_;      // Heat transfer coefficients
    mfem::Vector convectionTinf_;   // Ambient temperatures

    // Joule heating source
    std::unique_ptr<JouleHeatCoefficient> jouleHeatSource_;
    bool hasJouleHeating_ = false;

    std::unique_ptr<mfem::BilinearForm> aForm_;
    std::unique_ptr<mfem::LinearForm> bForm_;
    mfem::Array<int> essTdof_;

    const PhysicsProblemModel* problemModel_ = nullptr;

    void fillMaterialVectors(const PhysicsMaterialDatabase& materials);
    void buildConvectionBoundaryMarkers();
    void parseSourceTerms();
};

void HeatTransferSolver::Impl::fillMaterialVectors(const PhysicsMaterialDatabase& materials)
{
    const int maxAttribute = mesh_->attributes.Max();
    thermalK_.SetSize(maxAttribute);
    thermalK_ = 1.0;

    for (int attribute = 1; attribute <= maxAttribute; ++attribute) {
        auto tagIter = problemModel_->domainMaterialTag.find(attribute);
        if (tagIter == problemModel_->domainMaterialTag.end()) {
            continue;
        }

        auto propIter = materials.byTag.find(tagIter->second);
        if (propIter == materials.byTag.end()) {
            continue;
        }

        const MaterialPropertyModel& prop = propIter->second;
        thermalK_(attribute - 1) = prop.thermalConductivity > 0.0 
            ? prop.thermalConductivity : 1.0;
    }
}

void HeatTransferSolver::Impl::buildConvectionBoundaryMarkers()
{
    const int maxBdrAttr = mesh_->bdr_attributes.Max();
    convectionBdr_.SetSize(maxBdrAttr);
    convectionBdr_ = 0;
    convectionH_.SetSize(maxBdrAttr);
    convectionH_ = 0.0;
    convectionTinf_.SetSize(maxBdrAttr);
    convectionTinf_ = DEFAULT_REFERENCE_TEMPERATURE;

    for (const auto& bc : problemModel_->boundaries) {
        if (bc.field != FieldKind::Temperature) {
            continue;
        }
        if (bc.kind != BoundaryConditionKind::Convection) {
            continue;
        }

        for (int id : bc.boundaryIds) {
            if (id > 0 && id <= maxBdrAttr) {
                convectionBdr_[id - 1] = 1;
                convectionH_(id - 1) = bc.value;
                convectionTinf_(id - 1) = bc.auxiliary;
            }
        }
    }
}

void HeatTransferSolver::Impl::parseSourceTerms()
{
    hasJouleHeating_ = false;

    for (const auto& source : problemModel_->sources) {
        if (source.field != FieldKind::Temperature) {
            continue;
        }
        if (source.coupled && source.couplingKind == CouplingKind::JouleHeating) {
            hasJouleHeating_ = true;
            break;
        }
    }
}

HeatTransferSolver::HeatTransferSolver()
    : impl_(std::make_unique<Impl>())
{
}

HeatTransferSolver::~HeatTransferSolver() = default;

void HeatTransferSolver::setOrder(int order)
{
    impl_->order_ = order > 0 ? order : 1;
}

void HeatTransferSolver::setSolver(std::unique_ptr<LinearSolverStrategy> solver)
{
    impl_->solver_ = std::move(solver);
}

void HeatTransferSolver::initialize(mfem::Mesh& mesh,
                                    const PhysicsProblemModel& problemModel,
                                    const PhysicsMaterialDatabase& materials)
{
    impl_->mesh_ = &mesh;
    impl_->problemModel_ = &problemModel;

    // Create finite element collection and space
    impl_->fec_ = std::make_unique<mfem::H1_FECollection>(impl_->order_, mesh.Dimension());
    impl_->space_ = std::make_unique<mfem::FiniteElementSpace>(&mesh, impl_->fec_.get());

    // Initialize grid function
    impl_->temperature_ = std::make_unique<mfem::GridFunction>(impl_->space_.get());
    *impl_->temperature_ = DEFAULT_REFERENCE_TEMPERATURE;

    // Fill material vectors
    impl_->fillMaterialVectors(materials);

    // Create thermal conductivity coefficient
    impl_->thermalConductivity_ = std::make_unique<AttributeScalarCoefficient>();
    impl_->thermalConductivity_->setValues(impl_->thermalK_);

    // Build boundary markers
    impl_->buildConvectionBoundaryMarkers();

    // Parse source terms
    impl_->parseSourceTerms();

    // Create Joule heating source coefficient
    if (impl_->hasJouleHeating_) {
        impl_->jouleHeatSource_ = std::make_unique<JouleHeatCoefficient>();
    }
}

void HeatTransferSolver::applyBoundaryConditions()
{
    // Heat transfer has no essential boundary conditions in this model
    // (only convection on boundaries)
}

void HeatTransferSolver::assemble()
{
    impl_->aForm_ = std::make_unique<mfem::BilinearForm>(impl_->space_.get());
    impl_->aForm_->AddDomainIntegrator(
        new mfem::DiffusionIntegrator(*impl_->thermalConductivity_));

    // Add convection boundary condition: h * T term
    if (impl_->convectionBdr_.Max() > 0) {
        mfem::PWConstCoefficient hCoef(impl_->convectionH_);
        impl_->aForm_->AddBoundaryIntegrator(
            new mfem::MassIntegrator(hCoef), impl_->convectionBdr_);
    }

    impl_->aForm_->Assemble();

    // Build right-hand side
    impl_->bForm_ = std::make_unique<mfem::LinearForm>(impl_->space_.get());

    // Add Joule heating source
    if (impl_->hasJouleHeating_ && impl_->jouleHeatSource_) {
        impl_->bForm_->AddDomainIntegrator(
            new mfem::DomainLFIntegrator(*impl_->jouleHeatSource_));
    }

    // Add convection boundary term: h * T_inf
    if (impl_->convectionBdr_.Max() > 0) {
        mfem::PWConstCoefficient hCoef(impl_->convectionH_);
        mfem::PWConstCoefficient tinfCoef(impl_->convectionTinf_);
        mfem::ProductCoefficient htinfCoef(hCoef, tinfCoef);
        impl_->bForm_->AddBoundaryIntegrator(
            new mfem::BoundaryLFIntegrator(htinfCoef), impl_->convectionBdr_);
    }

    impl_->bForm_->Assemble();

    // No essential true DOFs for heat transfer in this model
    impl_->essTdof_.SetSize(0);
}

void HeatTransferSolver::solve()
{
    Check(impl_->solver_, "No linear solver set for HeatTransferSolver");

    mfem::SparseMatrix mat;
    mfem::Vector x, b;
    impl_->aForm_->FormLinearSystem(impl_->essTdof_,
                                    *impl_->temperature_,
                                    *impl_->bForm_,
                                    mat, x, b);
    x = 0.0;

    impl_->solver_->solve(mat, x, b);

    impl_->aForm_->RecoverFEMSolution(x, *impl_->bForm_, *impl_->temperature_);
}

mfem::GridFunction& HeatTransferSolver::getField()
{
    return *impl_->temperature_;
}

const mfem::GridFunction& HeatTransferSolver::getField() const
{
    return *impl_->temperature_;
}

FieldKind HeatTransferSolver::getFieldKind() const
{
    return FieldKind::Temperature;
}

mfem::FiniteElementSpace& HeatTransferSolver::getSpace()
{
    return *impl_->space_;
}

const mfem::FiniteElementSpace& HeatTransferSolver::getSpace() const
{
    return *impl_->space_;
}

void HeatTransferSolver::setPotentialField(const mfem::GridFunction* potential)
{
    if (impl_->jouleHeatSource_) {
        impl_->jouleHeatSource_->setPotential(potential);
    }
}

void HeatTransferSolver::setConductivityCoefficient(mfem::Coefficient* conductivity)
{
    if (impl_->jouleHeatSource_) {
        impl_->jouleHeatSource_->setConductivity(conductivity);
    }
}

mfem::Coefficient* HeatTransferSolver::getThermalConductivityCoefficient()
{
    return impl_->thermalConductivity_.get();
}

} // namespace mpfem
