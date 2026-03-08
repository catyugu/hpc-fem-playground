#ifndef MPFEM_PHYSICS_FIELD_SOLVER_HPP
#define MPFEM_PHYSICS_FIELD_SOLVER_HPP

#include "linear_solver_strategy.hpp"
#include "mpfem_types.hpp"
#include "physics_problem_model.hpp"

#include <memory>

namespace mpfem {

/**
 * @brief Abstract base class for physics field solvers.
 * 
 * Each physics field (electrostatics, heat transfer, solid mechanics) 
 * implements this interface, allowing unified coupling management.
 * 
 * Error handling: errors are logged and the program exits immediately.
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
     */
    virtual void initialize(FemMesh& mesh,
                           const PhysicsProblemModel& problemModel,
                           const PhysicsMaterialDatabase& materials) = 0;

    /**
     * @brief Apply boundary conditions to the system.
     */
    virtual void applyBoundaryConditions() = 0;

    /**
     * @brief Assemble the system matrix and right-hand side.
     */
    virtual void assemble() = 0;

    /**
     * @brief Solve the linear system.
     */
    virtual void solve() = 0;

    /**
     * @brief Get the field solution as a GridFunction.
     */
    virtual FemGridFunction& getField() = 0;
    virtual const FemGridFunction& getField() const = 0;

    /**
     * @brief Get the field kind this solver handles.
     */
    virtual FieldKind getFieldKind() const = 0;

    /**
     * @brief Get the finite element space for this field.
     */
    virtual FemFEspace& getSpace() = 0;
    virtual const FemFEspace& getSpace() const = 0;
};

} // namespace mpfem

#endif