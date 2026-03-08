#ifndef MPFEM_CASE_MODEL_HPP
#define MPFEM_CASE_MODEL_HPP

#include <set>
#include <string>
#include <vector>

namespace mpfem {

/**
 * @brief Scalar variable definition from case XML.
 */
struct VariableEntry {
    std::string name;
    std::string valueText;
    double siValue = 0.0;
};

/**
 * @brief Domain to material tag mapping rule.
 */
struct MaterialAssignment {
    std::set<int> domainIds;
    std::string materialTag;
};

/**
 * @brief Boundary condition definition for one physics block.
 */
struct BoundaryCondition {
    std::string kind;
    std::set<int> ids;
    std::string valueText;
    std::string auxText;
};

/**
 * @brief Volumetric or surface source definition.
 */
struct SourceDefinition {
    std::string kind;
    std::set<int> domainIds;
    std::string valueText;
};

/**
 * @brief Physics block extracted from case XML.
 */
struct PhysicsDefinition {
    std::string kind;
    std::vector<BoundaryCondition> boundaries;
    std::vector<SourceDefinition> sources;
};

/**
 * @brief Cross-physics coupling definition extracted from case XML.
 */
struct CoupledPhysicsDefinition {
    std::string name;
    std::string kind;
    std::vector<std::string> physicsKinds;
    std::set<int> domainIds;
};

/**
 * @brief Fully parsed case configuration.
 */
struct CaseDefinition {
    std::string caseName;
    std::string studyType;
    std::string meshPath;
    std::string materialsPath;
    std::string comsolResultPath;
    std::vector<VariableEntry> variables;
    std::vector<MaterialAssignment> materialAssignments;
    std::vector<PhysicsDefinition> physicsDefinitions;
    std::vector<CoupledPhysicsDefinition> coupledPhysicsDefinitions;
};

} // namespace mpfem

#endif
