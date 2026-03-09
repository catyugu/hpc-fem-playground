#ifndef MPFEM_PHYSICS_PROBLEM_MODEL_HPP
#define MPFEM_PHYSICS_PROBLEM_MODEL_HPP

#include "case_model.hpp"
#include "material_model.hpp"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace mpfem {

enum class FieldKind {
    ElectricPotential,
    Temperature,
    Displacement
};

enum class CouplingKind {
    JouleHeating,
    ThermalExpansion,
    Unknown
};

/**
 * @brief Boundary condition parameters (resolved to numeric values).
 * 
 * Stores the kind of boundary condition and its resolved parameter values.
 * Each boundary id maps to one BoundaryParams instance.
 */
struct BoundaryParams {
    std::string kind;
    std::map<std::string, double> values;
};

/// Boundary conditions organized by boundary id for each physics field
using BoundaryConditions = std::map<int, BoundaryParams>;

struct FieldSource {
    FieldKind field = FieldKind::Temperature;
    std::set<int> domainIds;
    bool coupled = false;
    CouplingKind couplingKind = CouplingKind::Unknown;
    double value = 0.0;
};

struct FieldConfig {
    int order = 1;
    SolverConfig solver;
};

struct PhysicsProblemModel {
    std::string caseName;
    std::string studyType;
    std::string meshPath;
    std::string comsolResultPath;
    std::map<std::string, double> variables;
    std::map<int, std::string> domainMaterialTag;
    std::map<FieldKind, BoundaryConditions> boundaries;  // field -> id -> params
    std::vector<FieldSource> sources;
    std::vector<CouplingKind> couplings;
    std::map<FieldKind, FieldConfig> fieldConfigs;
    CouplingConfig couplingConfig;
};

struct MaterialPropertyModel {
    std::string tag;
    std::string label;
    std::map<std::string, double> properties;
    double youngModulus = 0.0;
    double poissonRatio = 0.0;
    double rho0 = 0.0;
    double alpha = 0.0;
    double tref = 293.15;
    double electricConductivity = 0.0;
    double thermalConductivity = 0.0;
    double thermalExpansion = 0.0;
};

struct PhysicsMaterialDatabase {
    std::map<std::string, MaterialPropertyModel> byTag;
};

} // namespace mpfem

#endif
