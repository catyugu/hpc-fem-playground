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

enum class BoundaryConditionKind {
    Dirichlet,
    Insulation,
    Convection,
    Fixed,
    Free
};

enum class CouplingKind {
    JouleHeating,
    ThermalExpansion,
    Unknown
};

struct FieldBoundaryCondition {
    FieldKind field = FieldKind::ElectricPotential;
    BoundaryConditionKind kind = BoundaryConditionKind::Insulation;
    std::set<int> boundaryIds;
    double value = 0.0;
    double auxiliary = 0.0;
};

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
    std::vector<FieldBoundaryCondition> boundaries;
    std::vector<FieldSource> sources;
    std::vector<CouplingKind> couplings;
    std::map<FieldKind, FieldConfig> fieldConfigs;
    CouplingConfig couplingConfig;
};

struct MaterialPropertyModel {
    std::string tag;
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
