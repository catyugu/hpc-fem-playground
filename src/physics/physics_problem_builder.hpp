#ifndef MPFEM_PHYSICS_PROBLEM_BUILDER_HPP
#define MPFEM_PHYSICS_PROBLEM_BUILDER_HPP

#include "case_model.hpp"
#include "material_model.hpp"
#include "physics_problem_model.hpp"

namespace mpfem {

/**
 * @brief Builds PhysicsProblemModel from case and material definitions.
 */
class PhysicsProblemBuilder {
public:
    /**
     * @brief Builds the physics problem model.
     */
    static void build(const CaseDefinition &caseDefinition,
                      const PhysicsMaterialDatabase &materials,
                      PhysicsProblemModel &problemModel);
};

} // namespace mpfem

#endif