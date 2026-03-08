#ifndef MPFEM_MFEM_COUPLED_SOLVER_HPP
#define MPFEM_MFEM_COUPLED_SOLVER_HPP

#include "physics_problem_model.hpp"
#include "simulation_data.hpp"
#include "mpfem_types.hpp"

#include <map>
#include <set>
#include <string>

namespace mpfem {

/**
 * @brief MFEM-based solver for electro-thermal-mechanical coupled system.
 */
class MfemCoupledSolver {
public:
    /**
     * @brief Runs segregated coupling using MFEM weak-form assembly.
     * @param mesh Loaded MFEM mesh.
     * @param problemModel Generic problem model (includes coupling config).
     * @param materials Reduced material property model.
     * @param result Output sampled field results.
     */
    void solve(FemMesh &mesh,
               const PhysicsProblemModel &problemModel,
               const PhysicsMaterialDatabase &materials,
               CoupledFieldResult &result) const;
};

} // namespace mpfem

#endif