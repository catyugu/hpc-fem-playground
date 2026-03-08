#include "electrostatics_solver.hpp"
#include "linear_solver_strategy.hpp"

#include <cmath>

namespace mpfem {

namespace {

const double MIN_RESISTIVITY = 1e-16;
const double DEFAULT_REFERENCE_TEMPERATURE = 298;

/**
 * @brief Conductivity coefficient that supports temperature dependence.
 */
class ConductivityCoefficient : public mfem::Coefficient {
public:
    ConductivityCoefficient()
        : rho0_(), alpha_(), tref_(), sigma0_(), temperature_(nullptr) {}

    void setMaterialFields(const mfem::Vector& rho0,
                          const mfem::Vector& alpha,
                          const mfem::Vector& tref,
                          const mfem::Vector& sigma0)
    {
        rho0_ = rho0;
        alpha_ = alpha;
        tref_ = tref;
        sigma0_ = sigma0;
    }

    void setTemperatureField(const mfem::GridFunction* temperature)
    {
        temperature_ = temperature;
    }

    double Eval(mfem::ElementTransformation& transformation,
                const mfem::IntegrationPoint& integrationPoint) override
    {
        const int attribute = transformation.Attribute;
        if (attribute <= 0 || attribute > sigma0_.Size()) {
            return 1.0;
        }

        const double rho0 = rho0_(attribute - 1);
        if (rho0 > 0.0) {
            const double alpha = alpha_(attribute - 1);
            const double tref = tref_(attribute - 1);
            double temperatureValue = tref;
            if (temperature_ != nullptr) {
                temperatureValue = temperature_->GetValue(transformation, integrationPoint);
            }
            const double rho = rho0 * (1.0 + alpha * (temperatureValue - tref));
            if (std::abs(rho) > MIN_RESISTIVITY) {
                return 1.0 / rho;
            }
            return 1.0 / MIN_RESISTIVITY;
        }

        const double sigma = sigma0_(attribute - 1);
        return sigma > 0.0 ? sigma : 1.0;
    }

private:
    mfem::Vector rho0_;
    mfem::Vector alpha_;
    mfem::Vector tref_;
    mfem::Vector sigma0_;
    const mfem::GridFunction* temperature_;
};

} // namespace

class ElectrostaticsSolver::Impl {
public:
    int order_ = 1;
    mfem::Mesh* mesh_ = nullptr;
    std::unique_ptr<mfem::H1_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> space_;
    std::unique_ptr<mfem::GridFunction> potential_;
    std::unique_ptr<LinearSolverStrategy> solver_;
    
    std::unique_ptr<ConductivityCoefficient> conductivity_;
    mfem::Vector rho0_, alpha_, tref_, sigma0_;
    
    mfem::Array<int> essentialBdr_;
    mfem::Vector dirichletValues_;
    
    std::unique_ptr<mfem::BilinearForm> aForm_;
    std::unique_ptr<mfem::LinearForm> bForm_;
    mfem::Array<int> essTdof_;
    
    const PhysicsProblemModel* problemModel_ = nullptr;

    void fillMaterialVectors(const PhysicsMaterialDatabase& materials);
    void buildBoundaryMarkers();
};

void ElectrostaticsSolver::Impl::fillMaterialVectors(const PhysicsMaterialDatabase& materials)
{
    const int maxAttribute = mesh_->attributes.Max();
    rho0_.SetSize(maxAttribute);
    alpha_.SetSize(maxAttribute);
    tref_.SetSize(maxAttribute);
    sigma0_.SetSize(maxAttribute);

    rho0_ = 0.0;
    alpha_ = 0.0;
    tref_ = DEFAULT_REFERENCE_TEMPERATURE;
    sigma0_ = 1.0;

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
        rho0_(attribute - 1) = prop.rho0;
        alpha_(attribute - 1) = prop.alpha;
        tref_(attribute - 1) = prop.tref;
        sigma0_(attribute - 1) = prop.electricConductivity > 0.0 
            ? prop.electricConductivity : 1.0;
    }
}

void ElectrostaticsSolver::Impl::buildBoundaryMarkers()
{
    const int maxBdrAttr = mesh_->bdr_attributes.Max();
    essentialBdr_.SetSize(maxBdrAttr);
    essentialBdr_ = 0;
    dirichletValues_.SetSize(maxBdrAttr);
    dirichletValues_ = 0.0;

    for (const auto& bc : problemModel_->boundaries) {
        if (bc.field != FieldKind::ElectricPotential) {
            continue;
        }
        if (bc.kind != BoundaryConditionKind::Dirichlet) {
            continue;
        }

        for (int id : bc.boundaryIds) {
            if (id > 0 && id <= maxBdrAttr) {
                essentialBdr_[id - 1] = 1;
                dirichletValues_(id - 1) = bc.value;
            }
        }
    }
}

ElectrostaticsSolver::ElectrostaticsSolver()
    : impl_(std::make_unique<Impl>())
{
}

ElectrostaticsSolver::~ElectrostaticsSolver() = default;

void ElectrostaticsSolver::setOrder(int order)
{
    impl_->order_ = order > 0 ? order : 1;
}

void ElectrostaticsSolver::setSolver(std::unique_ptr<LinearSolverStrategy> solver)
{
    impl_->solver_ = std::move(solver);
}

bool ElectrostaticsSolver::initialize(mfem::Mesh& mesh,
                                      const PhysicsProblemModel& problemModel,
                                      const PhysicsMaterialDatabase& materials,
                                      std::string& errorMessage)
{
    errorMessage.clear();
    impl_->mesh_ = &mesh;
    impl_->problemModel_ = &problemModel;

    // Create finite element collection and space
    impl_->fec_ = std::make_unique<mfem::H1_FECollection>(impl_->order_, mesh.Dimension());
    impl_->space_ = std::make_unique<mfem::FiniteElementSpace>(&mesh, impl_->fec_.get());

    // Initialize grid function
    impl_->potential_ = std::make_unique<mfem::GridFunction>(impl_->space_.get());
    *impl_->potential_ = 0.0;

    // Fill material vectors
    impl_->fillMaterialVectors(materials);

    // Create conductivity coefficient
    impl_->conductivity_ = std::make_unique<ConductivityCoefficient>();
    impl_->conductivity_->setMaterialFields(
        impl_->rho0_, impl_->alpha_, impl_->tref_, impl_->sigma0_);

    // Build boundary markers
    impl_->buildBoundaryMarkers();

    return true;
}

void ElectrostaticsSolver::applyBoundaryConditions()
{
    mfem::PWConstCoefficient boundaryCoef(impl_->dirichletValues_);
    impl_->potential_->ProjectBdrCoefficient(boundaryCoef, impl_->essentialBdr_);
}

bool ElectrostaticsSolver::assemble(std::string& errorMessage)
{
    errorMessage.clear();

    impl_->aForm_ = std::make_unique<mfem::BilinearForm>(impl_->space_.get());
    impl_->aForm_->AddDomainIntegrator(
        new mfem::DiffusionIntegrator(*impl_->conductivity_));
    impl_->aForm_->Assemble();

    impl_->bForm_ = std::make_unique<mfem::LinearForm>(impl_->space_.get());
    impl_->bForm_->Assemble();

    impl_->space_->GetEssentialTrueDofs(impl_->essentialBdr_, impl_->essTdof_);

    return true;
}

bool ElectrostaticsSolver::solve(std::string& errorMessage)
{
    errorMessage.clear();

    mfem::SparseMatrix mat;
    mfem::Vector x, b;
    impl_->aForm_->FormLinearSystem(impl_->essTdof_, 
                                    *impl_->potential_, 
                                    *impl_->bForm_, 
                                    mat, x, b);

    if (!impl_->solver_) {
        errorMessage = "No linear solver set for ElectrostaticsSolver";
        return false;
    }

    if (!impl_->solver_->solve(mat, x, b, errorMessage)) {
        return false;
    }

    impl_->aForm_->RecoverFEMSolution(x, *impl_->bForm_, *impl_->potential_);
    return true;
}

mfem::GridFunction& ElectrostaticsSolver::getField()
{
    return *impl_->potential_;
}

const mfem::GridFunction& ElectrostaticsSolver::getField() const
{
    return *impl_->potential_;
}

FieldKind ElectrostaticsSolver::getFieldKind() const
{
    return FieldKind::ElectricPotential;
}

mfem::FiniteElementSpace& ElectrostaticsSolver::getSpace()
{
    return *impl_->space_;
}

const mfem::FiniteElementSpace& ElectrostaticsSolver::getSpace() const
{
    return *impl_->space_;
}

void ElectrostaticsSolver::setTemperatureField(const mfem::GridFunction* temperature)
{
    impl_->conductivity_->setTemperatureField(temperature);
}

mfem::Coefficient* ElectrostaticsSolver::getConductivityCoefficient()
{
    return impl_->conductivity_.get();
}

} // namespace mpfem
