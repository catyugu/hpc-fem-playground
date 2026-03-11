#ifndef MPFEM_HEAT_TRANSFER_SOLVER_HPP
#define MPFEM_HEAT_TRANSFER_SOLVER_HPP

#include "physics_field_solver.hpp"
#include "linear_solver_strategy.hpp"
#include "mfem.hpp"

#include <memory>

namespace mpfem {

/**
 * @brief Attribute-based scalar coefficient for material properties.
 */
class AttributeCoefficient : public mfem::Coefficient {
public:
    explicit AttributeCoefficient() = default;
    void setValues(const mfem::Vector& values) { values_ = values; }
    double Eval(mfem::ElementTransformation& tr, const mfem::IntegrationPoint& ip) override;

private:
    mfem::Vector values_;
};

/**
 * @brief Joule heating coefficient: Q = sigma * |grad(V)|^2
 */
class JouleHeatCoefficient : public mfem::Coefficient {
public:
    void setPotential(const mfem::GridFunction* potential) { potential_ = potential; }
    void setConductivity(mfem::Coefficient* conductivity) { conductivity_ = conductivity; }
    double Eval(mfem::ElementTransformation& tr, const mfem::IntegrationPoint& ip) override;

private:
    const mfem::GridFunction* potential_ = nullptr;
    mfem::Coefficient* conductivity_ = nullptr;
};

/**
 * @brief Solver for heat transfer (temperature) field.
 * 
 * Solves: -div(k * grad(T)) = Q
 * with support for convection BCs and Joule heating source.
 */
class HeatTransferSolver : public PhysicsFieldSolver {
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

    void setPotentialField(const mfem::GridFunction* potential);
    void setConductivityCoefficient(mfem::Coefficient* conductivity);
    mfem::Coefficient* getThermalConductivityCoefficient();

private:
    void fillMaterialVectors(const PhysicsMaterialDatabase& materials);
    void buildConvectionBoundaryMarkers();
    void parseSourceTerms();

    int order_ = 1;
    mfem::Mesh* mesh_ = nullptr;
    const PhysicsProblemModel* problemModel_ = nullptr;

    std::unique_ptr<mfem::H1_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> space_;
    std::unique_ptr<mfem::GridFunction> temperature_;
    std::unique_ptr<LinearSolverStrategy> solver_;

    // Thermal conductivity (attribute-based)
    AttributeCoefficient thermalConductivity_;
    mfem::Vector thermalK_;

    // Convection boundary conditions
    mfem::Array<int> convectionBdr_;
    mfem::Vector convectionH_;
    mfem::Vector convectionTinf_;

    // Convection coefficients (constant, created once)
    std::unique_ptr<mfem::PWConstCoefficient> hCoef_;
    std::unique_ptr<mfem::PWConstCoefficient> tinfCoef_;
    std::unique_ptr<mfem::ProductCoefficient> htinfCoef_;

    // Joule heating
    JouleHeatCoefficient jouleHeatSource_;
    bool hasJouleHeating_ = false;

    // Forms (created once, reassembled each iteration)
    std::unique_ptr<mfem::BilinearForm> aForm_;
    std::unique_ptr<mfem::LinearForm> bForm_;
    mfem::Array<int> essTdof_;  // Empty for heat transfer (no Dirichlet BCs)
};

} // namespace mpfem

#endif
