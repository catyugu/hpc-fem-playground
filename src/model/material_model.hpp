#ifndef MPFEM_MATERIAL_MODEL_HPP
#define MPFEM_MATERIAL_MODEL_HPP

#include <map>
#include <string>
#include <vector>

namespace mpfem {

/**
 * @brief Material with selected SI properties used by baseline solvers.
 */
struct MaterialDefinition {
    std::string tag;
    std::string label;
    std::map<std::string, double> siProperties;
};

/**
 * @brief Collection of materials loaded from material XML.
 */
struct MaterialDatabase {
    std::vector<MaterialDefinition> materials;
};

} // namespace mpfem

#endif
