#ifndef MPFEM_BOUNDARY_SPEC_HPP
#define MPFEM_BOUNDARY_SPEC_HPP

#include <map>
#include <string>
#include <vector>

namespace mpfem {

/**
 * @brief Specification for a boundary condition type.
 * 
 * Defines the required and optional parameters for a specific boundary condition kind.
 */
struct BoundarySpec {
    std::string kind;
    std::vector<std::string> requiredParams;
    std::vector<std::string> optionalParams;
};

/**
 * @brief Registry for boundary condition specifications.
 * 
 * Provides validation and parameter resolution for boundary conditions.
 * Each physics field can have its own set of boundary condition types.
 */
class BoundarySpecRegistry {
public:
    /// Register a boundary condition specification
    static void registerSpec(const std::string& physicsKind, const BoundarySpec& spec);
    
    /// Get specification for a boundary condition type
    static const BoundarySpec* getSpec(const std::string& physicsKind, const std::string& bcKind);
    
    /// Check if a boundary condition type is registered
    static bool hasSpec(const std::string& physicsKind, const std::string& bcKind);
    
    /// Initialize default specifications for all physics types
    static void initializeDefaults();

private:
    static std::map<std::string, std::map<std::string, BoundarySpec>> specs_;
};

} // namespace mpfem

#endif
