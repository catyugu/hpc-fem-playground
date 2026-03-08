#include "mfem_coupled_solver.hpp"

#include "coupling_manager.hpp"
#include "electrostatics_solver.hpp"
#include "heat_transfer_solver.hpp"
#include "linear_solver_strategy.hpp"
#include "solid_mechanics_solver.hpp"

#include <cmath>

namespace mpfem {

namespace {

bool recoverVertexSample(const mfem::Mesh& mesh,
                         const mfem::GridFunction& potential,
                         const mfem::GridFunction& temperature,
                         const mfem::GridFunction& displacement,
                         CoupledFieldResult& result,
                         std::string& errorMessage)
{
    errorMessage.clear();

    const int vertexCount = mesh.GetNV();
    result.coordinates.resize(vertexCount);
    result.electricPotential.resize(vertexCount, 0.0);
    result.temperature.resize(vertexCount, 0.0);
    result.displacement.resize(vertexCount, 0.0);

    mfem::Vector nodalPotential;
    mfem::Vector nodalTemperature;
    mfem::Vector nodalUx;
    mfem::Vector nodalUy;
    mfem::Vector nodalUz;

    potential.GetNodalValues(nodalPotential);
    temperature.GetNodalValues(nodalTemperature);
    displacement.GetNodalValues(nodalUx, 1);
    displacement.GetNodalValues(nodalUy, 2);
    displacement.GetNodalValues(nodalUz, 3);

    if (nodalPotential.Size() < vertexCount
        || nodalTemperature.Size() < vertexCount
        || nodalUx.Size() < vertexCount
        || nodalUy.Size() < vertexCount
        || nodalUz.Size() < vertexCount) {
        errorMessage = "Nodal sample size mismatch while exporting coupled result";
        return false;
    }

    for (int i = 0; i < vertexCount; ++i) {
        const double* vertex = mesh.GetVertex(i);
        result.coordinates[i].x = vertex[0];
        result.coordinates[i].y = vertex[1];
        result.coordinates[i].z = vertex[2];
        result.electricPotential[i] = nodalPotential(i);
        result.temperature[i] = nodalTemperature(i);
        const double ux = nodalUx(i);
        const double uy = nodalUy(i);
        const double uz = nodalUz(i);
        result.displacement[i] = std::sqrt(ux * ux + uy * uy + uz * uz);
    }

    return true;
}

FieldKind mapPhysicsKindToField(const std::string& physicsKind)
{
    if (physicsKind == "electrostatics") {
        return FieldKind::ElectricPotential;
    }
    if (physicsKind == "heat_transfer") {
        return FieldKind::Temperature;
    }
    return FieldKind::Displacement;
}

} // namespace

bool MfemCoupledSolver::solve(mfem::Mesh& mesh,
                              const PhysicsProblemModel& problemModel,
                              const PhysicsMaterialDatabase& materials,
                              int maxIterations,
                              double tolerance,
                              CoupledFieldResult& result,
                              std::string& errorMessage) const
{
    errorMessage.clear();
    result = CoupledFieldResult();

    // Create physics field solvers
    auto electrostaticsSolver = std::make_unique<ElectrostaticsSolver>();
    auto heatTransferSolver = std::make_unique<HeatTransferSolver>();
    auto solidMechanicsSolver = std::make_unique<SolidMechanicsSolver>();

    // Get field configurations from problem model
    auto getSolverConfig = [&](FieldKind kind) -> SolverConfig {
        auto iter = problemModel.fieldConfigs.find(kind);
        if (iter != problemModel.fieldConfigs.end()) {
            return iter->second.solver;
        }
        return SolverConfig{};
    };

    auto getFieldOrder = [&](FieldKind kind) -> int {
        auto iter = problemModel.fieldConfigs.find(kind);
        if (iter != problemModel.fieldConfigs.end()) {
            return iter->second.order;
        }
        return 1;
    };

    // Configure electrostatics solver
    int electroOrder = getFieldOrder(FieldKind::ElectricPotential);
    SolverConfig electroConfig = getSolverConfig(FieldKind::ElectricPotential);
    electrostaticsSolver->setOrder(electroOrder);
    electrostaticsSolver->setSolver(createLinearSolver(
        electroConfig.type, electroConfig.maxIterations, electroConfig.relativeTolerance));
    if (!electrostaticsSolver->initialize(mesh, problemModel, materials, errorMessage)) {
        return false;
    }

    // Configure heat transfer solver
    int thermalOrder = getFieldOrder(FieldKind::Temperature);
    SolverConfig thermalConfig = getSolverConfig(FieldKind::Temperature);
    heatTransferSolver->setOrder(thermalOrder);
    heatTransferSolver->setSolver(createLinearSolver(
        thermalConfig.type, thermalConfig.maxIterations, thermalConfig.relativeTolerance));
    if (!heatTransferSolver->initialize(mesh, problemModel, materials, errorMessage)) {
        return false;
    }

    // Configure solid mechanics solver
    int mechOrder = getFieldOrder(FieldKind::Displacement);
    SolverConfig mechConfig = getSolverConfig(FieldKind::Displacement);
    solidMechanicsSolver->setOrder(mechOrder);
    solidMechanicsSolver->setSolver(createLinearSolver(
        mechConfig.type, mechConfig.maxIterations, mechConfig.relativeTolerance));
    if (!solidMechanicsSolver->initialize(mesh, problemModel, materials, errorMessage)) {
        return false;
    }

    // Create coupling manager
    CouplingManager couplingManager;
    couplingManager.registerField(FieldKind::ElectricPotential, electrostaticsSolver.get());
    couplingManager.registerField(FieldKind::Temperature, heatTransferSolver.get());
    couplingManager.registerField(FieldKind::Displacement, solidMechanicsSolver.get());

    // Set coupling configuration from problem model, with fallback to function parameters
    CouplingConfig couplingConfig = problemModel.couplingConfig;
    if (maxIterations > 0) {
        couplingConfig.maxIterations = maxIterations;
    }
    if (tolerance > 0) {
        couplingConfig.tolerance = tolerance;
    }
    couplingManager.setCouplingConfig(couplingConfig);

    // Setup coupling
    couplingManager.setupCoupling();

    // Run coupled simulation
    if (!couplingManager.run(errorMessage)) {
        return false;
    }

    // Extract results
    const mfem::GridFunction& potential = electrostaticsSolver->getField();
    const mfem::GridFunction& temperature = heatTransferSolver->getField();
    const mfem::GridFunction& displacement = solidMechanicsSolver->getField();

    if (!recoverVertexSample(mesh, potential, temperature, displacement, result, errorMessage)) {
        return false;
    }

    return true;
}

} // namespace mpfem