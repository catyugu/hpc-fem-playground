#ifndef MPFEM_SOLID_MECHANICS_SOLVER_HPP
#define MPFEM_SOLID_MECHANICS_SOLVER_HPP

#include "physics_field_solver.hpp"

#include <memory>

namespace mpfem {

/**
 * @brief Solver for solid mechanics (displacement) field.
 * 
 * Solves: div(sigma) = 0  (static equilibrium)
 * with support for:
 * - Fixed constraint boundary conditions
 * - Thermal expansion load: sigma = (lambda + 2*mu) * alpha * deltaT * I
 */
class SolidMechanicsSolver : public PhysicsFieldSolver {
public:
    SolidMechanicsSolver();
    ~SolidMechanicsSolver() override;

    void setOrder(int order) override;
    void setSolver(std::unique_ptr<LinearSolverStrategy> solver) override;

    bool initialize(mfem::Mesh& mesh,
                   const PhysicsProblemModel& problemModel,
                   const PhysicsMaterialDatabase& materials,
                   std::string& errorMessage) override;

    void applyBoundaryConditions() override;
    bool assemble(std::string& errorMessage) override;
    bool solve(std::string& errorMessage) override;

    mfem::GridFunction& getField() override;
    const mfem::GridFunction& getField() const override;
    FieldKind getFieldKind() const override;
    mfem::FiniteElementSpace& getSpace() override;
    const mfem::FiniteElementSpace& getSpace() const override;

    /**
     * @brief Set the temperature field for thermal expansion load.
     * @param temperature Temperature field grid function.
     */
    void setTemperatureField(const mfem::GridFunction* temperature);

    /**
     * @brief Get the Lambda coefficient (Lame's first parameter).
     */
    mfem::Coefficient* getLambdaCoefficient();

    /**
     * @brief Get the Mu coefficient (Lame's second parameter / shear modulus).
     */
    mfem::Coefficient* getMuCoefficient();

    /**
     * @brief Get the thermal expansion coefficient.
     */
    mfem::Coefficient* getThermalExpansionCoefficient();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace mpfem

#endif
