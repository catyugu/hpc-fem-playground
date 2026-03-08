#include "physics_problem_builder.hpp"

#include "value_parser.hpp"

namespace mpfem {

namespace {

FieldKind mapPhysicsKindToField(const std::string &physicsKind)
{
    if (physicsKind == "electrostatics") {
        return FieldKind::ElectricPotential;
    }
    if (physicsKind == "heat_transfer") {
        return FieldKind::Temperature;
    }
    return FieldKind::Displacement;
}

BoundaryConditionKind mapBoundaryKind(const std::string &boundaryKind)
{
    if (boundaryKind == "voltage" || boundaryKind == "fixed_constraint") {
        return BoundaryConditionKind::Dirichlet;
    }
    if (boundaryKind == "electric_insulation" || boundaryKind == "thermal_insulation") {
        return BoundaryConditionKind::Insulation;
    }
    if (boundaryKind == "convection") {
        return BoundaryConditionKind::Convection;
    }
    if (boundaryKind == "free_boundary") {
        return BoundaryConditionKind::Free;
    }

    return BoundaryConditionKind::Insulation;
}

CouplingKind mapCouplingKind(const std::string &couplingKind)
{
    if (couplingKind == "joule_heating") {
        return CouplingKind::JouleHeating;
    }
    if (couplingKind == "thermal_expansion") {
        return CouplingKind::ThermalExpansion;
    }

    return CouplingKind::Unknown;
}

double getValueOrDefault(const std::map<std::string, double> &values,
                         const std::string &name,
                         const double defaultValue)
{
    const std::map<std::string, double>::const_iterator iter = values.find(name);
    if (iter == values.end()) {
        double parsedValue = 0.0;
        std::string errorMessage;
        if (ValueParser::parseFirstNumber(name, parsedValue, errorMessage)) {
            return parsedValue;
        }
        return defaultValue;
    }
    return iter->second;
}

} // namespace

bool PhysicsProblemBuilder::build(const CaseDefinition &caseDefinition,
                                  const MaterialDatabase &materialDatabase,
                                  PhysicsProblemModel &problemModel,
                                  PhysicsMaterialDatabase &physicsMaterials,
                                  std::string &errorMessage)
{
    errorMessage.clear();
    problemModel = PhysicsProblemModel();
    physicsMaterials = PhysicsMaterialDatabase();

    problemModel.caseName = caseDefinition.caseName;
    problemModel.studyType = caseDefinition.studyType;
    problemModel.meshPath = caseDefinition.meshPath;
    problemModel.comsolResultPath = caseDefinition.comsolResultPath;

    for (std::size_t i = 0; i < caseDefinition.variables.size(); ++i) {
        problemModel.variables[caseDefinition.variables[i].name] = caseDefinition.variables[i].siValue;
    }

    for (std::size_t i = 0; i < caseDefinition.materialAssignments.size(); ++i) {
        const MaterialAssignment &assignment = caseDefinition.materialAssignments[i];
        for (std::set<int>::const_iterator domainIter = assignment.domainIds.begin();
             domainIter != assignment.domainIds.end();
             ++domainIter) {
            problemModel.domainMaterialTag[*domainIter] = assignment.materialTag;
        }
    }

    for (std::size_t i = 0; i < caseDefinition.physicsDefinitions.size(); ++i) {
        const PhysicsDefinition &physicsDefinition = caseDefinition.physicsDefinitions[i];
        const FieldKind fieldKind = mapPhysicsKindToField(physicsDefinition.kind);

        for (std::size_t j = 0; j < physicsDefinition.boundaries.size(); ++j) {
            const BoundaryCondition &boundary = physicsDefinition.boundaries[j];
            FieldBoundaryCondition modelBoundary;
            modelBoundary.field = fieldKind;
            modelBoundary.kind = mapBoundaryKind(boundary.kind);
            modelBoundary.boundaryIds = boundary.ids;
            modelBoundary.value = getValueOrDefault(problemModel.variables, boundary.valueText, 0.0);
            modelBoundary.auxiliary = getValueOrDefault(problemModel.variables, boundary.auxText, 0.0);
            problemModel.boundaries.push_back(modelBoundary);
        }

        for (std::size_t j = 0; j < physicsDefinition.sources.size(); ++j) {
            const SourceDefinition &source = physicsDefinition.sources[j];
            FieldSource modelSource;
            modelSource.field = fieldKind;
            modelSource.domainIds = source.domainIds;
            modelSource.value = getValueOrDefault(problemModel.variables, source.valueText, 0.0);
            modelSource.coupled = source.valueText == "joule_heating";
            modelSource.couplingKind = modelSource.coupled ? CouplingKind::JouleHeating : CouplingKind::Unknown;
            problemModel.sources.push_back(modelSource);
        }
    }

    for (std::size_t i = 0; i < caseDefinition.coupledPhysicsDefinitions.size(); ++i) {
        problemModel.couplings.push_back(mapCouplingKind(caseDefinition.coupledPhysicsDefinitions[i].kind));
    }

    for (std::size_t i = 0; i < materialDatabase.materials.size(); ++i) {
        const MaterialDefinition &material = materialDatabase.materials[i];
        MaterialPropertyModel property;
        property.tag = material.tag;
        property.youngModulus = getValueOrDefault(material.siProperties, "E", 1.0);
        property.poissonRatio = getValueOrDefault(material.siProperties, "nu", 0.3);
        property.rho0 = getValueOrDefault(material.siProperties, "rho0", 0.0);
        property.alpha = getValueOrDefault(material.siProperties, "alpha", 0.0);
        property.tref = getValueOrDefault(material.siProperties, "Tref", 293.15);
        property.electricConductivity = getValueOrDefault(material.siProperties, "electricconductivity", 1.0);
        property.thermalConductivity = getValueOrDefault(material.siProperties, "thermalconductivity", 1.0);
        property.thermalExpansion = getValueOrDefault(material.siProperties, "thermalexpansioncoefficient", 1.0e-5);
        physicsMaterials.byTag[property.tag] = property;
    }

    if (problemModel.domainMaterialTag.empty()) {
        errorMessage = "No domain-material mapping found while building physics model";
        return false;
    }

    return true;
}

} // namespace mpfem
