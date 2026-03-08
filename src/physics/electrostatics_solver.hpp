#ifndef MPFEM_ELECTROSTATICS_SOLVER_HPP
#define MPFEM_ELECTROSTATICS_SOLVER_HPP

#include "physics_field_solver.hpp"

#include <memory>

namespace mpfem {

/**
 * @brief Solver for electrostatics (electric potential) field.
 * 
 * Solves: -div(sigma * grad(V)) = 0
 * where sigma is the electrical conductivity (can be temperature-dependent).
 */
class ElectrostaticsSolver : public PhysicsFieldSolver {
public:
    ElectrostaticsSolver();
    ~ElectrostaticsSolver() override;

    void setOrder(int order) override;
    void setSolver(std::unique_ptr<LinearSolverStrategy> solver) override;

    void initialize(FemMesh& mesh,
                   const PhysicsProblemModel& problemModel,
                   const PhysicsMaterialDatabase& materials) override;

    void applyBoundaryConditions() override;
    void assemble() override;
    void solve() override;

    FemGridFunction& getField() override;
    const FemGridFunction& getField() const override;
    FieldKind getFieldKind() const override;
    FemFEspace& getSpace() override;
    const FemFEspace& getSpace() const override;

    /**
     * @brief Set the temperature field for temperature-dependent conductivity.
     */
    void setTemperatureField(const FemGridFunction* temperature);

    /**
     * @brief Get the conductivity coefficient (for coupling).
     */
    mfem::Coefficient* getConductivityCoefficient();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace mpfem

#endif