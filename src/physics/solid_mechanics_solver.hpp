#ifndef MPFEM_SOLID_MECHANICS_SOLVER_HPP
#define MPFEM_SOLID_MECHANICS_SOLVER_HPP

#include "physics_field_solver.hpp"
#include "heat_transfer_solver.hpp"  // For AttributeCoefficient
#include "linear_solver_strategy.hpp"
#include "mfem.hpp"

#include <memory>
#include <optional>

namespace mpfem {

/**
 * @brief Thermal expansion gradient load coefficient.
 */
class ThermalExpansionLoadCoefficient : public mfem::VectorCoefficient {
public:
    explicit ThermalExpansionLoadCoefficient(int dimension)
        : mfem::VectorCoefficient(dimension * dimension)
        , dimension_(dimension)
    {}

    void setTemperatureField(const mfem::GridFunction* temperature) { temperature_ = temperature; }
    void setAlphaCoefficient(mfem::Coefficient* alpha) { alpha_ = alpha; }
    void setLambdaCoefficient(mfem::Coefficient* lambda) { lambda_ = lambda; }
    void setMuCoefficient(mfem::Coefficient* mu) { mu_ = mu; }
    void setReferenceTemperature(double tref) { referenceTemperature_ = tref; }

    void Eval(mfem::Vector& V, mfem::ElementTransformation& tr,
              const mfem::IntegrationPoint& ip) override;

private:
    int dimension_;
    const mfem::GridFunction* temperature_ = nullptr;
    mfem::Coefficient* alpha_ = nullptr;
    mfem::Coefficient* lambda_ = nullptr;
    mfem::Coefficient* mu_ = nullptr;
    double referenceTemperature_ = 293.15;
};

/**
 * @brief Solver for solid mechanics (displacement) field.
 * 
 * Solves linear elasticity with thermal expansion.
 */
class SolidMechanicsSolver : public PhysicsFieldSolver {
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
    mfem::Coefficient* getLambdaCoefficient();
    mfem::Coefficient* getMuCoefficient();
    mfem::Coefficient* getThermalExpansionCoefficient();

private:
    void fillMaterialVectors(const PhysicsMaterialDatabase& materials);
    void computeLameParameters();
    void buildBoundaryMarkers();
    void parseCouplingTerms();

    int order_ = 1;
    mfem::Mesh* mesh_ = nullptr;
    const PhysicsProblemModel* problemModel_ = nullptr;

    std::unique_ptr<mfem::H1_FECollection> fec_;
    std::unique_ptr<mfem::FiniteElementSpace> space_;
    std::unique_ptr<mfem::GridFunction> displacement_;
    std::unique_ptr<LinearSolverStrategy> solver_;

    // Material coefficients
    mfem::Vector young_, poisson_, thermalExpansion_;
    AttributeCoefficient lambda_;
    AttributeCoefficient mu_;
    AttributeCoefficient alpha_;
    mfem::Vector lambdaValues_, muValues_;

    // Thermal expansion
    std::optional<ThermalExpansionLoadCoefficient> thermalExpansionLoad_;
    bool hasThermalExpansion_ = false;
    double referenceTemperature_ = 293.15;

    // Boundary conditions
    mfem::Array<int> fixedBdr_;
    mfem::Array<int> essTdof_;

    // Zero displacement coefficient for fixed boundaries
    mfem::Vector zeroDisp_;
    std::unique_ptr<mfem::VectorConstantCoefficient> zeroCoef_;

    // Forms
    std::unique_ptr<mfem::BilinearForm> aForm_;
    std::unique_ptr<mfem::LinearForm> bForm_;
};

} // namespace mpfem

#endif
