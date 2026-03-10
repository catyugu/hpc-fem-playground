#include "physics_problem_builder.hpp"
#include "boundary_spec.hpp"
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

double resolveValue(const std::map<std::string, double> &variables,
                    const std::string &text,
                    double defaultValue = 0.0)
{
    if (text.empty()) {
        return defaultValue;
    }
    // First try to find in variables map
    auto iter = variables.find(text);
    if (iter != variables.end()) {
        return iter->second;
    }
    // Otherwise parse as number
    double result = defaultValue;
    ValueParser::parseFirstNumber(text, result);
    return result;
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

        // Initialize boundary conditions map for this field
        BoundaryConditions& bcMap = problemModel.boundaries[fieldKind];

        // Process each boundary condition definition
        for (const auto& boundary : physicsDef.boundaries) {
            // Resolve parameters to numeric values
            std::map<std::string, double> resolvedParams;
            for (const auto& [name, text] : boundary.params) {
                resolvedParams[name] = resolveValue(problemModel.variables, text);
            }

            // Store boundary condition for each boundary id
            for (int id : boundary.ids) {
                BoundaryParams params;
                params.kind = boundary.kind;
                params.values = resolvedParams;
                bcMap[id] = params;
            }
        }

        for (const auto& source : physicsDef.sources) {
            FieldSource modelSource;
            modelSource.field = fieldKind;
            modelSource.domainIds = source.domainIds;
            modelSource.value = resolveValue(problemModel.variables, source.valueText);
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
