#include "physics_problem_builder.hpp"
#include "value_parser.hpp"
#include "logger.hpp"

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
                         double defaultValue)
{
    auto iter = values.find(name);
    if (iter == values.end()) {
        double result = defaultValue;
        ValueParser::parseFirstNumber(name, result);
        return result;
    }
    return iter->second;
}

} // namespace

void PhysicsProblemBuilder::build(const CaseDefinition &caseDefinition,
                                  const PhysicsMaterialDatabase &materials,
                                  PhysicsProblemModel &problemModel)
{
    problemModel = PhysicsProblemModel();

    problemModel.caseName = caseDefinition.caseName;
    problemModel.studyType = caseDefinition.studyType;
    problemModel.meshPath = caseDefinition.meshPath;
    problemModel.comsolResultPath = caseDefinition.comsolResultPath;

    for (const auto& variable : caseDefinition.variables) {
        problemModel.variables[variable.name] = variable.siValue;
    }

    for (const auto& assignment : caseDefinition.materialAssignments) {
        for (int domain : assignment.domainIds) {
            problemModel.domainMaterialTag[domain] = assignment.materialTag;
        }
    }

    for (const auto& physicsDef : caseDefinition.physicsDefinitions) {
        const FieldKind fieldKind = mapPhysicsKindToField(physicsDef.kind);

        // Store field configuration
        FieldConfig fieldConfig;
        fieldConfig.order = physicsDef.order;
        fieldConfig.solver = physicsDef.solver;
        problemModel.fieldConfigs[fieldKind] = fieldConfig;

        for (const auto& boundary : physicsDef.boundaries) {
            FieldBoundaryCondition modelBoundary;
            modelBoundary.field = fieldKind;
            modelBoundary.kind = mapBoundaryKind(boundary.kind);
            modelBoundary.boundaryIds = boundary.ids;
            modelBoundary.value = getValueOrDefault(problemModel.variables, boundary.valueText, 0.0);
            modelBoundary.auxiliary = getValueOrDefault(problemModel.variables, boundary.auxText, 0.0);
            problemModel.boundaries.push_back(modelBoundary);
        }

        for (const auto& source : physicsDef.sources) {
            FieldSource modelSource;
            modelSource.field = fieldKind;
            modelSource.domainIds = source.domainIds;
            modelSource.value = getValueOrDefault(problemModel.variables, source.valueText, 0.0);
            modelSource.coupled = source.valueText == "joule_heating";
            modelSource.couplingKind = modelSource.coupled ? CouplingKind::JouleHeating : CouplingKind::Unknown;
            problemModel.sources.push_back(modelSource);
        }
    }

    for (const auto& couplingDef : caseDefinition.coupledPhysicsDefinitions) {
        problemModel.couplings.push_back(mapCouplingKind(couplingDef.kind));
    }

    problemModel.couplingConfig = caseDefinition.couplingConfig;

    Check(!problemModel.domainMaterialTag.empty(), 
          "No domain-material mapping found while building physics model");
}

} // namespace mpfem