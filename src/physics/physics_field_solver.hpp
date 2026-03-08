#ifndef MPFEM_PHYSICS_FIELD_SOLVER_HPP
#define MPFEM_PHYSICS_FIELD_SOLVER_HPP

#include "linear_solver_strategy.hpp"
#include "physics_problem_model.hpp"

#include "mfem.hpp"

#include <memory>
#include <string>

namespace mpfem {

/**
 * @brief Abstract base class for physics field solvers.
 * 
 * Each physics field (electrostatics, heat transfer, solid mechanics) 
 * implements this interface, allowing unified coupling management.
 */
class PhysicsFieldSolver {
public:
    virtual ~PhysicsFieldSolver() = default;

    /**
     * @brief Set the polynomial order for this field.
     */
    virtual void setOrder(int order) = 0;

    /**
     * @brief Set the linear solver strategy.
     */
    virtual void setSolver(std::unique_ptr<LinearSolverStrategy> solver) = 0;

    /**
     * @brief Initialize the solver with mesh and materials.
     * @param mesh The finite element mesh.
     * @param problemModel The physics problem model containing boundary conditions.
     * @param materials The material database.
     * @param errorMessage Error message on failure.
     * @return True if initialization succeeded.
     */
    virtual bool initialize(mfem::Mesh& mesh,
                           const PhysicsProblemModel& problemModel,
                           const PhysicsMaterialDatabase& materials,
                           std::string& errorMessage) = 0;

    /**
     * @brief Apply boundary conditions to the system.
     */
    virtual void applyBoundaryConditions() = 0;

    /**
     * @brief Assemble the system matrix and right-hand side.
     * @param errorMessage Error message on failure.
     * @return True if assembly succeeded.
     */
    virtual bool assemble(std::string& errorMessage) = 0;

    /**
     * @brief Solve the linear system.
     * @param errorMessage Error message on failure.
     * @return True if solve succeeded.
     */
    virtual bool solve(std::string& errorMessage) = 0;

    /**
     * @brief Get the field solution as a GridFunction.
     */
    virtual mfem::GridFunction& getField() = 0;
    virtual const mfem::GridFunction& getField() const = 0;

    /**
     * @brief Get the field kind this solver handles.
     */
    virtual FieldKind getFieldKind() const = 0;

    /**
     * @brief Get the finite element space for this field.
     */
    virtual mfem::FiniteElementSpace& getSpace() = 0;
    virtual const mfem::FiniteElementSpace& getSpace() const = 0;
};

} // namespace mpfem

#endif
