#ifndef MPFEM_SOLID_MECHANICS_SOLVER_HPP
#define MPFEM_SOLID_MECHANICS_SOLVER_HPP

#include "physics_field_solver.hpp"

#include <memory>

namespace mpfem {

/**
 * @brief Solver for solid mechanics (displacement) field.
 * 
 * Solves linear elasticity with thermal expansion:
 * -div(sigma) = 0
 * where sigma = lambda * div(u) * I + 2 * mu * epsilon(u) - (3*lambda + 2*mu) * alpha * deltaT * I
 */
class SolidMechanicsSolver : public PhysicsFieldSolver {
public:
    SolidMechanicsSolver();
    ~SolidMechanicsSolver() override;

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
     * @brief Set the temperature field for thermal expansion.
     */
    void setTemperatureField(const mfem::GridFunction* temperature);

    /**
     * @brief Get Lame parameter lambda coefficient.
     */
    mfem::Coefficient* getLambdaCoefficient();

    /**
     * @brief Get Lame parameter mu coefficient.
     */
    mfem::Coefficient* getMuCoefficient();

    /**
     * @brief Get thermal expansion coefficient.
     */
    mfem::Coefficient* getThermalExpansionCoefficient();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace mpfem

#endif