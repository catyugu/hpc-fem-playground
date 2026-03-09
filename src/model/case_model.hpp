#ifndef MPFEM_CASE_MODEL_HPP
#define MPFEM_CASE_MODEL_HPP

#include <set>
#include <string>
#include <vector>

namespace mpfem {

/**
 * @brief Coupling iteration method.
 */
enum class CouplingMethod {
    NewtonRaphson,
    Picard
};

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
 * @brief Linear solver configuration.
 */
struct SolverConfig {
    std::string type = "cg_gs";  // "cg_gs", "gmres"
    int maxIterations = 1000;
    double relativeTolerance = 1e-10;
    int printLevel = 0;
};

/**
 * @brief Physics block extracted from case XML.
 */
struct PhysicsDefinition {
    std::string kind;
    int order = 1;  // Polynomial order for this physics field
    SolverConfig solver;
    std::vector<BoundaryCondition> boundaries;
    std::vector<SourceDefinition> sources;
};

/**
 * @brief Coupling iteration configuration.
 */
struct CouplingConfig {
    std::string method = "newton_raphson";  // "newton_raphson", "picard"
    int maxIterations = 15;
    double tolerance = 1e-6;
};

/**
 * @brief Convert method string to enum.
 */
inline CouplingMethod parseCouplingMethod(const std::string& method) {
    if (method == "picard") return CouplingMethod::Picard;
    return CouplingMethod::NewtonRaphson;
}

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
    CouplingConfig couplingConfig;
};

} // namespace mpfem

#endif
