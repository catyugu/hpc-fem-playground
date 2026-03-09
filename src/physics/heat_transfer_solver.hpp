#ifndef MPFEM_HEAT_TRANSFER_SOLVER_HPP
#define MPFEM_HEAT_TRANSFER_SOLVER_HPP

#include "physics_field_solver.hpp"

#include <memory>

namespace mpfem {

/**
 * @brief Solver for heat transfer (temperature) field.
 * 
 * Solves: -div(k * grad(T)) = Q
 * with support for:
 * - Convection boundary conditions: k * grad(T) . n = h * (T - T_inf)
 * - Joule heating source: Q = sigma * |grad(V)|^2
 */
class HeatTransferSolver : public PhysicsFieldSolver {
public:
    HeatTransferSolver();
    ~HeatTransferSolver() override;

    void setOrder(int order) override;
    void setSolver(std::unique_ptr<LinearSolverStrategy> solver) override;

    void initialize(mfem::Mesh& mesh,
                   const PhysicsProblemModel& problemModel,
                   const PhysicsMaterialDatabase& materials) override;

    void applyBoundaryConditions() override;
    void assemble() override;
    void solve() override;

    mfem::GridFunction& getField() override;
    const mfem::GridFunction& getField() const override;
    FieldKind getFieldKind() const override;
    mfem::FiniteElementSpace& getSpace() override;
    const mfem::FiniteElementSpace& getSpace() const override;

    /**
     * @brief Set the electric potential field for Joule heating source.
     */
    void setPotentialField(const mfem::GridFunction* potential);

    /**
     * @brief Set the conductivity coefficient for Joule heating calculation.
     */
    void setConductivityCoefficient(mfem::Coefficient* conductivity);

    /**
     * @brief Get the thermal conductivity coefficient (for coupling).
     */
    mfem::Coefficient* getThermalConductivityCoefficient();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace mpfem

#endif