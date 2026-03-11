#ifndef MPFEM_ELECTROSTATICS_SOLVER_HPP
#define MPFEM_ELECTROSTATICS_SOLVER_HPP

#include "physics_field_solver.hpp"
#include "linear_solver_strategy.hpp"
#include "mfem.hpp"

#include <memory>

namespace mpfem {

/**
 * @brief Coefficient for temperature-dependent electrical conductivity.
 */
class ConductivityCoefficient : public mfem::Coefficient {
public:
    void setMaterialFields(const mfem::Vector& rho0,
                           const mfem::Vector& alpha,
                           const mfem::Vector& tref,
                           const mfem::Vector& sigma0);
    void setTemperatureField(const mfem::GridFunction* temperature);
    double Eval(mfem::ElementTransformation& transformation,
                const mfem::IntegrationPoint& ip) override;

private:
    mfem::Vector rho0_, alpha_, tref_, sigma0_;
    const mfem::GridFunction* temperature_ = nullptr;
};

/**
 * @brief Solver for electrostatics (electric potential) field.
 * 
 * Solves: -div(sigma * grad(V)) = 0
 * where sigma is the electrical conductivity (can be temperature-dependent).
 */
class ElectrostaticsSolver : public PhysicsFieldSolver {
public:
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

    void setTemperatureField(const mfem::GridFunction* temperature);
    mfem::Coefficient* getConductivityCoefficient();

private:
    void fillMaterialVectors(const PhysicsMaterialDatabase& materials);
    void buildBoundaryMarkers();

    int order_ = 1;
    mfem::Mesh* mesh_ = nullptr;
    const PhysicsProblemModel* problemModel_ = nullptr;

    std::unique_ptr<mfem::H1_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> space_;
    std::unique_ptr<mfem::GridFunction> potential_;
    std::unique_ptr<LinearSolverStrategy> solver_;

    // Material coefficients
    ConductivityCoefficient conductivity_;
    mfem::Vector rho0_, alpha_, tref_, sigma0_;

    // Boundary conditions
    mfem::Array<int> essentialBdr_;
    mfem::Vector dirichletValues_;
    mfem::Array<int> essTdof_;  // Cached essential true DOFs

    // Forms (created once in initialize, reassembled each iteration)
    std::unique_ptr<mfem::BilinearForm> aForm_;
    std::unique_ptr<mfem::LinearForm> bForm_;
};

} // namespace mpfem

#endif
