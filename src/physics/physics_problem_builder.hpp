#ifndef MPFEM_PHYSICS_PROBLEM_BUILDER_HPP
#define MPFEM_PHYSICS_PROBLEM_BUILDER_HPP

#include "case_model.hpp"
#include "material_model.hpp"
#include "physics_problem_model.hpp"

#include <string>

namespace mpfem {

/**
 * @brief Builds configuration-independent physics problem models.
 */
class PhysicsProblemBuilder {
public:
    /**
     * @brief Converts parsed case and material data to generic physics model.
     * @param caseDefinition Parsed case definition.
     * @param materialDatabase Parsed material database.
     * @param problemModel Output generic problem model.
     * @param physicsMaterials Output reduced material property model.
     * @param errorMessage Error details on failure.
     * @return True when conversion succeeds.
     */
    static bool build(const CaseDefinition &caseDefinition,
                      const MaterialDatabase &materialDatabase,
                      PhysicsProblemModel &problemModel,
                      PhysicsMaterialDatabase &physicsMaterials,
                      std::string &errorMessage);
};

} // namespace mpfem

#endif
