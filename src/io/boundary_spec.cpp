#include "boundary_spec.hpp"

namespace mpfem {

std::map<std::string, std::map<std::string, BoundarySpec>> BoundarySpecRegistry::specs_;

void BoundarySpecRegistry::registerSpec(const std::string& physicsKind, const BoundarySpec& spec)
{
    specs_[physicsKind][spec.kind] = spec;
}

const BoundarySpec* BoundarySpecRegistry::getSpec(const std::string& physicsKind, const std::string& bcKind)
{
    auto physicsIt = specs_.find(physicsKind);
    if (physicsIt == specs_.end()) {
        return nullptr;
    }
    auto bcIt = physicsIt->second.find(bcKind);
    if (bcIt == physicsIt->second.end()) {
        return nullptr;
    }
    return &bcIt->second;
}

bool BoundarySpecRegistry::hasSpec(const std::string& physicsKind, const std::string& bcKind)
{
    return getSpec(physicsKind, bcKind) != nullptr;
}

void BoundarySpecRegistry::initializeDefaults()
{
    // Electrostatics boundary conditions
    registerSpec("electrostatics", {"voltage", {"value"}, {}});
    registerSpec("electrostatics", {"electric_insulation", {}, {}});
    registerSpec("electrostatics", {"ground", {"value"}, {}});
    
    // Heat transfer boundary conditions
    registerSpec("heat_transfer", {"temperature", {"value"}, {}});
    registerSpec("heat_transfer", {"thermal_insulation", {}, {}});
    registerSpec("heat_transfer", {"convection", {"h", "T_inf"}, {}});
    registerSpec("heat_transfer", {"heat_flux", {"q"}, {}});
    
    // Solid mechanics boundary conditions
    registerSpec("solid_mechanics", {"fixed_constraint", {}, {}});
    registerSpec("solid_mechanics", {"displacement", {"value"}, {}});
    registerSpec("solid_mechanics", {"free_boundary", {}, {}});
    registerSpec("solid_mechanics", {"surface_traction", {"tx", "ty", "tz"}, {"tz"}});
    registerSpec("solid_mechanics", {"pressure", {"p"}, {}});
}

} // namespace mpfem
