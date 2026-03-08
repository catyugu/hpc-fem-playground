#include "mfem_coupled_solver.hpp"

#include "linear_system_solver.hpp"

#include <cmath>

namespace mpfem {

namespace {

const double MIN_RESISTIVITY = 1e-16;
const double DEFAULT_REFERENCE_TEMPERATURE = 293.15;
const int FIELD_ORDER = 1;

class AttributeScalarCoefficient : public mfem::Coefficient {
public:
    AttributeScalarCoefficient() : values_() {}

    void setValues(const mfem::Vector &values)
    {
        values_ = values;
    }

    double Eval(mfem::ElementTransformation &transformation,
                const mfem::IntegrationPoint &integrationPoint) override
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

class ConductivityCoefficient : public mfem::Coefficient {
public:
    ConductivityCoefficient()
        : rho0_(), alpha_(), tref_(), sigma0_(), temperature_(NULL) {}

    void setMaterialFields(const mfem::Vector &rho0,
                           const mfem::Vector &alpha,
                           const mfem::Vector &tref,
                           const mfem::Vector &sigma0)
    {
        rho0_ = rho0;
        alpha_ = alpha;
        tref_ = tref;
        sigma0_ = sigma0;
    }

    void setTemperatureField(const mfem::GridFunction *temperature)
    {
        temperature_ = temperature;
    }

    double Eval(mfem::ElementTransformation &transformation,
                const mfem::IntegrationPoint &integrationPoint) override
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
            if (temperature_ != NULL) {
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
    const mfem::GridFunction *temperature_;
};

class JouleHeatCoefficient : public mfem::Coefficient {
public:
    JouleHeatCoefficient()
        : potential_(NULL), conductivity_(NULL), gradient_() {}

    void setPotential(const mfem::GridFunction *potential)
    {
        potential_ = potential;
    }

    void setConductivity(mfem::Coefficient *conductivity)
    {
        conductivity_ = conductivity;
    }

    double Eval(mfem::ElementTransformation &transformation,
                const mfem::IntegrationPoint &integrationPoint) override
    {
        if (potential_ == NULL || conductivity_ == NULL) {
            return 0.0;
        }

        transformation.SetIntPoint(&integrationPoint);
        potential_->GetGradient(transformation, gradient_);
        const double sigma = conductivity_->Eval(transformation, integrationPoint);
        const double norm2 = gradient_ * gradient_;
        return sigma * norm2;
    }

private:
    const mfem::GridFunction *potential_;
    mfem::Coefficient *conductivity_;
    mutable mfem::Vector gradient_;
};

class ThermalExpansionGradLoadCoefficient : public mfem::VectorCoefficient {
public:
        ThermalExpansionGradLoadCoefficient(const int dimension)
                : mfem::VectorCoefficient(dimension * dimension),
                    dimension_(dimension),
          temperature_(NULL),
          alpha_(NULL),
                    lambda_(NULL),
                    mu_(NULL),
                    referenceTemperature_(DEFAULT_REFERENCE_TEMPERATURE),
                    loadVector_() {}

    void setTemperatureField(const mfem::GridFunction *temperature)
    {
        temperature_ = temperature;
    }

    void setAlphaCoefficient(mfem::Coefficient *alpha)
    {
        alpha_ = alpha;
    }

    void setLambdaCoefficient(mfem::Coefficient *lambda)
    {
        lambda_ = lambda;
    }

    void setMuCoefficient(mfem::Coefficient *mu)
    {
        mu_ = mu;
    }

    void setReferenceTemperature(const double referenceTemperature)
    {
        referenceTemperature_ = referenceTemperature;
    }

    void Eval(mfem::Vector &vector,
              mfem::ElementTransformation &transformation,
              const mfem::IntegrationPoint &integrationPoint) override
    {
        vector.SetSize(dimension_ * dimension_);
        vector = 0.0;
        if (temperature_ == NULL || alpha_ == NULL || lambda_ == NULL || mu_ == NULL) {
            return;
        }

        transformation.SetIntPoint(&integrationPoint);
        const double temperatureValue = temperature_->GetValue(transformation, integrationPoint);
        const double alphaValue = alpha_->Eval(transformation, integrationPoint);
        const double lambdaValue = lambda_->Eval(transformation, integrationPoint);
        const double muValue = mu_->Eval(transformation, integrationPoint);
        const double deltaTemperature = temperatureValue - referenceTemperature_;

        const double thermalModulus = (3.0 * lambdaValue + 2.0 * muValue) * alphaValue * deltaTemperature;

        loadVector_.SetSize(dimension_ * dimension_);
        loadVector_ = 0.0;
        for (int i = 0; i < dimension_; ++i) {
            loadVector_(i * dimension_ + i) = thermalModulus;
        }
        vector = loadVector_;
    }

private:
    int dimension_;
    const mfem::GridFunction *temperature_;
    mfem::Coefficient *alpha_;
    mfem::Coefficient *lambda_;
    mfem::Coefficient *mu_;
    double referenceTemperature_;
    mutable mfem::Vector loadVector_;
};

void fillMaterialVectors(const mfem::Mesh &mesh,
                         const PhysicsProblemModel &problemModel,
                         const PhysicsMaterialDatabase &materials,
                         mfem::Vector &rho0,
                         mfem::Vector &alpha,
                         mfem::Vector &tref,
                         mfem::Vector &sigma0,
                         mfem::Vector &thermalK,
                         mfem::Vector &young,
                         mfem::Vector &poisson,
                         mfem::Vector &thermalExpansion)
{
    const int maxAttribute = mesh.attributes.Max();
    rho0.SetSize(maxAttribute);
    alpha.SetSize(maxAttribute);
    tref.SetSize(maxAttribute);
    sigma0.SetSize(maxAttribute);
    thermalK.SetSize(maxAttribute);
    young.SetSize(maxAttribute);
    poisson.SetSize(maxAttribute);
    thermalExpansion.SetSize(maxAttribute);

    rho0 = 0.0;
    alpha = 0.0;
    tref = DEFAULT_REFERENCE_TEMPERATURE;
    sigma0 = 1.0;
    thermalK = 1.0;
    young = 1.0;
    poisson = 0.3;
    thermalExpansion = 1.0e-5;

    for (int attribute = 1; attribute <= maxAttribute; ++attribute) {
        const std::map<int, std::string>::const_iterator materialTagIter =
            problemModel.domainMaterialTag.find(attribute);
        if (materialTagIter == problemModel.domainMaterialTag.end()) {
            continue;
        }

        const std::map<std::string, MaterialPropertyModel>::const_iterator propertyIter =
            materials.byTag.find(materialTagIter->second);
        if (propertyIter == materials.byTag.end()) {
            continue;
        }

        const MaterialPropertyModel &property = propertyIter->second;
        rho0(attribute - 1) = property.rho0;
        alpha(attribute - 1) = property.alpha;
        tref(attribute - 1) = property.tref;
        sigma0(attribute - 1) = property.electricConductivity > 0.0 ? property.electricConductivity : 1.0;
        thermalK(attribute - 1) = property.thermalConductivity > 0.0 ? property.thermalConductivity : 1.0;
        young(attribute - 1) = property.youngModulus > 0.0 ? property.youngModulus : 1.0;
        poisson(attribute - 1) = property.poissonRatio;
        thermalExpansion(attribute - 1) = property.thermalExpansion > 0.0 ? property.thermalExpansion : 1.0e-5;
    }
}

void buildBoundaryMarkers(const PhysicsProblemModel &problemModel,
                          const mfem::Mesh &mesh,
                          const FieldKind field,
                          const BoundaryConditionKind kind,
                          mfem::Array<int> &marker,
                          mfem::Vector &values,
                          const bool withValues)
{
    const int maxBoundaryAttribute = mesh.bdr_attributes.Max();
    marker.SetSize(maxBoundaryAttribute);
    marker = 0;
    if (withValues) {
        values.SetSize(maxBoundaryAttribute);
        values = 0.0;
    }

    for (std::size_t i = 0; i < problemModel.boundaries.size(); ++i) {
        const FieldBoundaryCondition &condition = problemModel.boundaries[i];
        if (condition.field != field || condition.kind != kind) {
            continue;
        }

        for (std::set<int>::const_iterator idIter = condition.boundaryIds.begin();
             idIter != condition.boundaryIds.end();
             ++idIter) {
            const int boundaryAttribute = *idIter;
            if (boundaryAttribute <= 0 || boundaryAttribute > maxBoundaryAttribute) {
                continue;
            }
            marker[boundaryAttribute - 1] = 1;
            if (withValues) {
                values(boundaryAttribute - 1) = condition.value;
            }
        }
    }
}

bool recoverVertexSample(const mfem::Mesh &mesh,
                         const mfem::GridFunction &potential,
                         const mfem::GridFunction &temperature,
                         const mfem::GridFunction &displacement,
                         CoupledFieldResult &result,
                         std::string &errorMessage)
{
    errorMessage.clear();

    const int vertexCount = mesh.GetNV();
    result.coordinates.resize(vertexCount);
    result.electricPotential.resize(vertexCount, 0.0);
    result.temperature.resize(vertexCount, 0.0);
    result.displacement.resize(vertexCount, 0.0);

    mfem::Vector nodalPotential;
    mfem::Vector nodalTemperature;
    mfem::Vector nodalUx;
    mfem::Vector nodalUy;
    mfem::Vector nodalUz;

    potential.GetNodalValues(nodalPotential);
    temperature.GetNodalValues(nodalTemperature);
    displacement.GetNodalValues(nodalUx, 1);
    displacement.GetNodalValues(nodalUy, 2);
    displacement.GetNodalValues(nodalUz, 3);

    if (nodalPotential.Size() < vertexCount
        || nodalTemperature.Size() < vertexCount
        || nodalUx.Size() < vertexCount
        || nodalUy.Size() < vertexCount
        || nodalUz.Size() < vertexCount) {
        errorMessage = "Nodal sample size mismatch while exporting coupled result";
        return false;
    }

    for (int i = 0; i < vertexCount; ++i) {
        const double *vertex = mesh.GetVertex(i);
        result.coordinates[i].x = vertex[0];
        result.coordinates[i].y = vertex[1];
        result.coordinates[i].z = vertex[2];
        result.electricPotential[i] = nodalPotential(i);
        result.temperature[i] = nodalTemperature(i);
        const double ux = nodalUx(i);
        const double uy = nodalUy(i);
        const double uz = nodalUz(i);
        result.displacement[i] = std::sqrt(ux * ux + uy * uy + uz * uz);
    }

    return true;
}

double computeFieldDelta(const mfem::GridFunction &previous,
                         const mfem::GridFunction &current)
{
    if (previous.Size() != current.Size()) {
        return 1.0;
    }

    double maxDelta = 0.0;
    for (int i = 0; i < previous.Size(); ++i) {
        const double delta = std::abs(previous(i) - current(i));
        if (delta > maxDelta) {
            maxDelta = delta;
        }
    }

    return maxDelta;
}

} // namespace

bool MfemCoupledSolver::solve(mfem::Mesh &mesh,
                              const PhysicsProblemModel &problemModel,
                              const PhysicsMaterialDatabase &materials,
                              int maxIterations,
                              double tolerance,
                              CoupledFieldResult &result,
                              std::string &errorMessage) const
{
    errorMessage.clear();
    result = CoupledFieldResult();

    mfem::H1_FECollection scalarCollection(FIELD_ORDER, mesh.Dimension());
    mfem::FiniteElementSpace scalarSpace(&mesh, &scalarCollection);

    mfem::H1_FECollection vectorCollection(FIELD_ORDER, mesh.Dimension());
    mfem::FiniteElementSpace vectorSpace(&mesh, &vectorCollection, mesh.Dimension());

    mfem::GridFunction potential(&scalarSpace);
    mfem::GridFunction temperature(&scalarSpace);
    mfem::GridFunction displacement(&vectorSpace);

    potential = 0.0;
    temperature = DEFAULT_REFERENCE_TEMPERATURE;
    displacement = 0.0;

    mfem::Vector rho0;
    mfem::Vector alpha;
    mfem::Vector tref;
    mfem::Vector sigma0;
    mfem::Vector thermalK;
    mfem::Vector young;
    mfem::Vector poisson;
    mfem::Vector thermalExpansion;

    fillMaterialVectors(mesh,
                        problemModel,
                        materials,
                        rho0,
                        alpha,
                        tref,
                        sigma0,
                        thermalK,
                        young,
                        poisson,
                        thermalExpansion);

    ConductivityCoefficient conductivityCoefficient;
    conductivityCoefficient.setMaterialFields(rho0, alpha, tref, sigma0);

    AttributeScalarCoefficient thermalConductivityCoefficient;
    thermalConductivityCoefficient.setValues(thermalK);

    mfem::Vector lambdaValues(young.Size());
    mfem::Vector muValues(young.Size());
    for (int i = 0; i < young.Size(); ++i) {
        const double e = young(i);
        const double nu = poisson(i);
        muValues(i) = e / (2.0 * (1.0 + nu));
        lambdaValues(i) = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    }

    AttributeScalarCoefficient lambdaCoefficient;
    lambdaCoefficient.setValues(lambdaValues);
    AttributeScalarCoefficient muCoefficient;
    muCoefficient.setValues(muValues);
    AttributeScalarCoefficient thermalExpansionCoefficient;
    thermalExpansionCoefficient.setValues(thermalExpansion);
    JouleHeatCoefficient jouleHeatCoefficient;
    jouleHeatCoefficient.setConductivity(&conductivityCoefficient);

    ThermalExpansionGradLoadCoefficient thermalExpansionLoadCoefficient(mesh.Dimension());
    thermalExpansionLoadCoefficient.setAlphaCoefficient(&thermalExpansionCoefficient);
    thermalExpansionLoadCoefficient.setLambdaCoefficient(&lambdaCoefficient);
    thermalExpansionLoadCoefficient.setMuCoefficient(&muCoefficient);
    thermalExpansionLoadCoefficient.setReferenceTemperature(DEFAULT_REFERENCE_TEMPERATURE);

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        mfem::GridFunction previousTemperature(temperature);
        mfem::GridFunction previousDisplacement(displacement);

        conductivityCoefficient.setTemperatureField(&temperature);

        mfem::Array<int> electricEssBdr;
        mfem::Vector electricDirichletValues;
        buildBoundaryMarkers(problemModel,
                             mesh,
                             FieldKind::ElectricPotential,
                             BoundaryConditionKind::Dirichlet,
                             electricEssBdr,
                             electricDirichletValues,
                             true);

        mfem::PWConstCoefficient electricBoundaryCoefficient(electricDirichletValues);
        potential.ProjectBdrCoefficient(electricBoundaryCoefficient, electricEssBdr);

        mfem::BilinearForm electroA(&scalarSpace);
        electroA.AddDomainIntegrator(new mfem::DiffusionIntegrator(conductivityCoefficient));
        mfem::LinearForm electroB(&scalarSpace);
        electroA.Assemble();
        electroB.Assemble();

        mfem::Array<int> electroEssTdof;
        scalarSpace.GetEssentialTrueDofs(electricEssBdr, electroEssTdof);
        mfem::SparseMatrix electroMat;
        mfem::Vector electroX;
        mfem::Vector electroRhs;
        electroA.FormLinearSystem(electroEssTdof, potential, electroB, electroMat, electroX, electroRhs);
        if (!LinearSystemSolver::solve(electroMat, electroRhs, electroX, 4000, 1e-10, errorMessage)) {
            return false;
        }
        electroA.RecoverFEMSolution(electroX, electroB, potential);

        jouleHeatCoefficient.setPotential(&potential);

        mfem::Array<int> convectionBdr;
        mfem::Vector convectionValues;
        buildBoundaryMarkers(problemModel,
                             mesh,
                             FieldKind::Temperature,
                             BoundaryConditionKind::Convection,
                             convectionBdr,
                             convectionValues,
                             true);

        mfem::Vector convectionAux(convectionValues.Size());
        convectionAux = 293.15;
        for (std::size_t i = 0; i < problemModel.boundaries.size(); ++i) {
            const FieldBoundaryCondition &condition = problemModel.boundaries[i];
            if (condition.field != FieldKind::Temperature
                || condition.kind != BoundaryConditionKind::Convection) {
                continue;
            }
            for (std::set<int>::const_iterator idIter = condition.boundaryIds.begin();
                 idIter != condition.boundaryIds.end();
                 ++idIter) {
                if (*idIter > 0 && *idIter <= convectionAux.Size()) {
                    convectionAux(*idIter - 1) = condition.auxiliary;
                }
            }
        }

        mfem::PWConstCoefficient hCoefficient(convectionValues);
        mfem::PWConstCoefficient tInfCoefficient(convectionAux);
        mfem::ProductCoefficient hTimesTinfCoefficient(hCoefficient, tInfCoefficient);

        mfem::BilinearForm thermalA(&scalarSpace);
        thermalA.AddDomainIntegrator(new mfem::DiffusionIntegrator(thermalConductivityCoefficient));
        thermalA.AddBoundaryIntegrator(new mfem::MassIntegrator(hCoefficient), convectionBdr);

        mfem::LinearForm thermalB(&scalarSpace);
        thermalB.AddDomainIntegrator(new mfem::DomainLFIntegrator(jouleHeatCoefficient));
        thermalB.AddBoundaryIntegrator(new mfem::BoundaryLFIntegrator(hTimesTinfCoefficient), convectionBdr);

        thermalA.Assemble();
        thermalB.Assemble();

        mfem::Array<int> thermalEssTdof;
        mfem::SparseMatrix thermalMat;
        mfem::Vector thermalX;
        mfem::Vector thermalRhs;
        thermalA.FormLinearSystem(thermalEssTdof, temperature, thermalB, thermalMat, thermalX, thermalRhs);
        if (!LinearSystemSolver::solve(thermalMat, thermalRhs, thermalX, 4000, 1e-10, errorMessage)) {
            return false;
        }
        thermalA.RecoverFEMSolution(thermalX, thermalB, temperature);

        thermalExpansionLoadCoefficient.setTemperatureField(&temperature);

        mfem::Array<int> fixedBdr;
        mfem::Vector dummyValues;
        buildBoundaryMarkers(problemModel,
                             mesh,
                             FieldKind::Displacement,
                             BoundaryConditionKind::Dirichlet,
                             fixedBdr,
                             dummyValues,
                             false);

        mfem::Vector zeroDisplacement(mesh.Dimension());
        zeroDisplacement = 0.0;
        mfem::VectorConstantCoefficient zeroCoefficient(zeroDisplacement);
        displacement.ProjectBdrCoefficient(zeroCoefficient, fixedBdr);

        mfem::BilinearForm mechanicsA(&vectorSpace);
        mechanicsA.AddDomainIntegrator(new mfem::ElasticityIntegrator(lambdaCoefficient, muCoefficient));

        mfem::LinearForm mechanicsB(&vectorSpace);
        mechanicsB.AddDomainIntegrator(new mfem::VectorDomainLFGradIntegrator(thermalExpansionLoadCoefficient));

        mechanicsA.Assemble();
        mechanicsB.Assemble();

        mfem::Array<int> mechanicsEssTdof;
        vectorSpace.GetEssentialTrueDofs(fixedBdr, mechanicsEssTdof);
        mfem::SparseMatrix mechanicsMat;
        mfem::Vector mechanicsX;
        mfem::Vector mechanicsRhs;
        mechanicsA.FormLinearSystem(mechanicsEssTdof,
                                    displacement,
                                    mechanicsB,
                                    mechanicsMat,
                                    mechanicsX,
                                    mechanicsRhs);
        if (!LinearSystemSolver::solve(mechanicsMat, mechanicsRhs, mechanicsX, 6000, 1e-9, errorMessage)) {
            return false;
        }
        mechanicsA.RecoverFEMSolution(mechanicsX, mechanicsB, displacement);

        const double temperatureDelta = computeFieldDelta(previousTemperature, temperature);
        const double displacementDelta = computeFieldDelta(previousDisplacement, displacement);
        if (temperatureDelta < tolerance && displacementDelta < tolerance) {
            break;
        }
    }

    if (!recoverVertexSample(mesh, potential, temperature, displacement, result, errorMessage)) {
        return false;
    }

    return true;
}

} // namespace mpfem
