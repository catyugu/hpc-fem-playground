#include "solid_mechanics_solver.hpp"
#include "linear_solver_strategy.hpp"

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
 * @brief Thermal expansion gradient load coefficient.
 * 
 * Computes the stress due to thermal expansion:
 * sigma_thermal = (3*lambda + 2*mu) * alpha * deltaT * I
 * which contributes as a body force: div(sigma_thermal)
 */
class ThermalExpansionGradLoadCoefficient : public mfem::VectorCoefficient {
public:
    explicit ThermalExpansionGradLoadCoefficient(int dimension)
        : mfem::VectorCoefficient(dimension * dimension)
        , dimension_(dimension)
        , temperature_(nullptr)
        , alpha_(nullptr)
        , lambda_(nullptr)
        , mu_(nullptr)
        , referenceTemperature_(DEFAULT_REFERENCE_TEMPERATURE)
        , loadVector_()
    {}

    void setTemperatureField(const mfem::GridFunction* temperature)
    {
        temperature_ = temperature;
    }

    void setAlphaCoefficient(mfem::Coefficient* alpha)
    {
        alpha_ = alpha;
    }

    void setLambdaCoefficient(mfem::Coefficient* lambda)
    {
        lambda_ = lambda;
    }

    void setMuCoefficient(mfem::Coefficient* mu)
    {
        mu_ = mu;
    }

    void setReferenceTemperature(double referenceTemperature)
    {
        referenceTemperature_ = referenceTemperature;
    }

    void Eval(mfem::Vector& vector,
              mfem::ElementTransformation& transformation,
              const mfem::IntegrationPoint& integrationPoint) override
    {
        vector.SetSize(dimension_ * dimension_);
        vector = 0.0;

        if (temperature_ == nullptr || alpha_ == nullptr 
            || lambda_ == nullptr || mu_ == nullptr) {
            return;
        }

        transformation.SetIntPoint(&integrationPoint);
        const double temperatureValue = temperature_->GetValue(transformation, integrationPoint);
        const double alphaValue = alpha_->Eval(transformation, integrationPoint);
        const double lambdaValue = lambda_->Eval(transformation, integrationPoint);
        const double muValue = mu_->Eval(transformation, integrationPoint);
        const double deltaTemperature = temperatureValue - referenceTemperature_;

        const double thermalModulus = (3.0 * lambdaValue + 2.0 * muValue) 
                                     * alphaValue * deltaTemperature;

        loadVector_.SetSize(dimension_ * dimension_);
        loadVector_ = 0.0;
        for (int i = 0; i < dimension_; ++i) {
            loadVector_(i * dimension_ + i) = thermalModulus;
        }
        vector = loadVector_;
    }

private:
    int dimension_;
    const mfem::GridFunction* temperature_;
    mfem::Coefficient* alpha_;
    mfem::Coefficient* lambda_;
    mfem::Coefficient* mu_;
    double referenceTemperature_;
    mutable mfem::Vector loadVector_;
};

} // namespace

class SolidMechanicsSolver::Impl {
public:
    int order_ = 1;
    mfem::Mesh* mesh_ = nullptr;
    std::unique_ptr<mfem::H1_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> space_;
    std::unique_ptr<mfem::GridFunction> displacement_;
    std::unique_ptr<LinearSolverStrategy> solver_;

    // Material coefficients
    mfem::Vector young_, poisson_, thermalExpansion_;
    std::unique_ptr<AttributeScalarCoefficient> lambda_;
    std::unique_ptr<AttributeScalarCoefficient> mu_;
    std::unique_ptr<AttributeScalarCoefficient> alpha_;
    mfem::Vector lambdaValues_, muValues_;

    // Thermal expansion load
    std::unique_ptr<ThermalExpansionGradLoadCoefficient> thermalExpansionLoad_;
    bool hasThermalExpansion_ = false;
    double referenceTemperature_ = DEFAULT_REFERENCE_TEMPERATURE;

    // Boundary conditions
    mfem::Array<int> fixedBdr_;

    std::unique_ptr<mfem::BilinearForm> aForm_;
    std::unique_ptr<mfem::LinearForm> bForm_;
    mfem::Array<int> essTdof_;

    const PhysicsProblemModel* problemModel_ = nullptr;

    void fillMaterialVectors(const PhysicsMaterialDatabase& materials);
    void computeLameParameters();
    void buildBoundaryMarkers();
    void parseCouplingTerms();
};

void SolidMechanicsSolver::Impl::fillMaterialVectors(const PhysicsMaterialDatabase& materials)
{
    const int maxAttribute = mesh_->attributes.Max();
    young_.SetSize(maxAttribute);
    poisson_.SetSize(maxAttribute);
    thermalExpansion_.SetSize(maxAttribute);

    young_ = 1.0;
    poisson_ = 0.3;
    thermalExpansion_ = 1.0e-5;

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
        young_(attribute - 1) = prop.youngModulus > 0.0 ? prop.youngModulus : 1.0;
        poisson_(attribute - 1) = prop.poissonRatio;
        thermalExpansion_(attribute - 1) = prop.thermalExpansion > 0.0 
            ? prop.thermalExpansion : 1.0e-5;
    }
}

void SolidMechanicsSolver::Impl::computeLameParameters()
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

    lambda_ = std::make_unique<AttributeScalarCoefficient>();
    lambda_->setValues(lambdaValues_);
    mu_ = std::make_unique<AttributeScalarCoefficient>();
    mu_->setValues(muValues_);
    alpha_ = std::make_unique<AttributeScalarCoefficient>();
    alpha_->setValues(thermalExpansion_);
}

void SolidMechanicsSolver::Impl::buildBoundaryMarkers()
{
    const int maxBdrAttr = mesh_->bdr_attributes.Max();
    fixedBdr_.SetSize(maxBdrAttr);
    fixedBdr_ = 0;

    for (const auto& bc : problemModel_->boundaries) {
        if (bc.field != FieldKind::Displacement) {
            continue;
        }
        if (bc.kind != BoundaryConditionKind::Dirichlet) {
            continue;
        }

        for (int id : bc.boundaryIds) {
            if (id > 0 && id <= maxBdrAttr) {
                fixedBdr_[id - 1] = 1;
            }
        }
    }
}

void SolidMechanicsSolver::Impl::parseCouplingTerms()
{
    hasThermalExpansion_ = false;

    for (const auto& coupling : problemModel_->couplings) {
        if (coupling == CouplingKind::ThermalExpansion) {
            hasThermalExpansion_ = true;
            break;
        }
    }
}

SolidMechanicsSolver::SolidMechanicsSolver()
    : impl_(std::make_unique<Impl>())
{
}

SolidMechanicsSolver::~SolidMechanicsSolver() = default;

void SolidMechanicsSolver::setOrder(int order)
{
    impl_->order_ = order > 0 ? order : 1;
}

void SolidMechanicsSolver::setSolver(std::unique_ptr<LinearSolverStrategy> solver)
{
    impl_->solver_ = std::move(solver);
}

bool SolidMechanicsSolver::initialize(mfem::Mesh& mesh,
                                      const PhysicsProblemModel& problemModel,
                                      const PhysicsMaterialDatabase& materials,
                                      std::string& errorMessage)
{
    errorMessage.clear();
    impl_->mesh_ = &mesh;
    impl_->problemModel_ = &problemModel;

    // Create vector finite element space for displacement
    impl_->fec_ = std::make_unique<mfem::H1_FECollection>(impl_->order_, mesh.Dimension());
    impl_->space_ = std::make_unique<mfem::FiniteElementSpace>(
        &mesh, impl_->fec_.get(), mesh.Dimension());

    // Initialize grid function
    impl_->displacement_ = std::make_unique<mfem::GridFunction>(impl_->space_.get());
    *impl_->displacement_ = 0.0;

    // Fill material vectors and compute Lame parameters
    impl_->fillMaterialVectors(materials);
    impl_->computeLameParameters();

    // Build boundary markers
    impl_->buildBoundaryMarkers();

    // Parse coupling terms
    impl_->parseCouplingTerms();

    // Create thermal expansion load coefficient if needed
    if (impl_->hasThermalExpansion_) {
        impl_->thermalExpansionLoad_ = std::make_unique<ThermalExpansionGradLoadCoefficient>(
            mesh.Dimension());
        impl_->thermalExpansionLoad_->setAlphaCoefficient(impl_->alpha_.get());
        impl_->thermalExpansionLoad_->setLambdaCoefficient(impl_->lambda_.get());
        impl_->thermalExpansionLoad_->setMuCoefficient(impl_->mu_.get());
        impl_->thermalExpansionLoad_->setReferenceTemperature(impl_->referenceTemperature_);
    }

    return true;
}

void SolidMechanicsSolver::applyBoundaryConditions()
{
    // Project zero displacement on fixed boundaries
    mfem::Vector zeroDisp(impl_->mesh_->Dimension());
    zeroDisp = 0.0;
    mfem::VectorConstantCoefficient zeroCoef(zeroDisp);
    impl_->displacement_->ProjectBdrCoefficient(zeroCoef, impl_->fixedBdr_);
}

bool SolidMechanicsSolver::assemble(std::string& errorMessage)
{
    errorMessage.clear();

    impl_->aForm_ = std::make_unique<mfem::BilinearForm>(impl_->space_.get());
    impl_->aForm_->AddDomainIntegrator(
        new mfem::ElasticityIntegrator(*impl_->lambda_, *impl_->mu_));
    impl_->aForm_->Assemble();

    // Build right-hand side
    impl_->bForm_ = std::make_unique<mfem::LinearForm>(impl_->space_.get());

    // Add thermal expansion load
    if (impl_->hasThermalExpansion_ && impl_->thermalExpansionLoad_) {
        impl_->bForm_->AddDomainIntegrator(
            new mfem::VectorDomainLFGradIntegrator(*impl_->thermalExpansionLoad_));
    }

    impl_->bForm_->Assemble();

    // Get essential true DOFs
    impl_->space_->GetEssentialTrueDofs(impl_->fixedBdr_, impl_->essTdof_);

    return true;
}

bool SolidMechanicsSolver::solve(std::string& errorMessage)
{
    errorMessage.clear();

    mfem::SparseMatrix mat;
    mfem::Vector x, b;
    impl_->aForm_->FormLinearSystem(impl_->essTdof_,
                                    *impl_->displacement_,
                                    *impl_->bForm_,
                                    mat, x, b);

    if (!impl_->solver_) {
        errorMessage = "No linear solver set for SolidMechanicsSolver";
        return false;
    }

    if (!impl_->solver_->solve(mat, x, b, errorMessage)) {
        return false;
    }

    impl_->aForm_->RecoverFEMSolution(x, *impl_->bForm_, *impl_->displacement_);
    return true;
}

mfem::GridFunction& SolidMechanicsSolver::getField()
{
    return *impl_->displacement_;
}

const mfem::GridFunction& SolidMechanicsSolver::getField() const
{
    return *impl_->displacement_;
}

FieldKind SolidMechanicsSolver::getFieldKind() const
{
    return FieldKind::Displacement;
}

mfem::FiniteElementSpace& SolidMechanicsSolver::getSpace()
{
    return *impl_->space_;
}

const mfem::FiniteElementSpace& SolidMechanicsSolver::getSpace() const
{
    return *impl_->space_;
}

void SolidMechanicsSolver::setTemperatureField(const mfem::GridFunction* temperature)
{
    if (impl_->thermalExpansionLoad_) {
        impl_->thermalExpansionLoad_->setTemperatureField(temperature);
    }
}

mfem::Coefficient* SolidMechanicsSolver::getLambdaCoefficient()
{
    return impl_->lambda_.get();
}

mfem::Coefficient* SolidMechanicsSolver::getMuCoefficient()
{
    return impl_->mu_.get();
}

mfem::Coefficient* SolidMechanicsSolver::getThermalExpansionCoefficient()
{
    return impl_->alpha_.get();
}

} // namespace mpfem
