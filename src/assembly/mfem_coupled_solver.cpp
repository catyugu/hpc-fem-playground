#include "mfem_coupled_solver.hpp"

#include "coupling_manager.hpp"
#include "electrostatics_solver.hpp"
#include "heat_transfer_solver.hpp"
#include "linear_solver_strategy.hpp"
#include "solid_mechanics_solver.hpp"
#include "logger.hpp"

#include <cmath>
#include <memory>
#include <vector>

namespace mpfem {

namespace {

void recoverVertexSample(const mfem::Mesh& mesh,
                         const mfem::GridFunction& potential,
                         const mfem::GridFunction& temperature,
                         const mfem::GridFunction& displacement,
                         CoupledFieldResult& result)
{
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

    Check(nodalPotential.Size() >= vertexCount, "Nodal potential size mismatch");
    Check(nodalTemperature.Size() >= vertexCount, "Nodal temperature size mismatch");
    Check(nodalUx.Size() >= vertexCount, "Nodal Ux size mismatch");
    Check(nodalUy.Size() >= vertexCount, "Nodal Uy size mismatch");
    Check(nodalUz.Size() >= vertexCount, "Nodal Uz size mismatch");

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
}

} // namespace

void MfemCoupledSolver::solve(mfem::Mesh& mesh,
                              const PhysicsProblemModel& problemModel,
                              const PhysicsMaterialDatabase& materials,
                              CoupledFieldResult& result) const
{
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
    electrostaticsSolver->initialize(mesh, problemModel, materials);

    // Configure heat transfer solver
    int thermalOrder = getFieldOrder(FieldKind::Temperature);
    SolverConfig thermalConfig = getSolverConfig(FieldKind::Temperature);
    heatTransferSolver->setOrder(thermalOrder);
    heatTransferSolver->setSolver(createLinearSolver(
        thermalConfig.type, thermalConfig.maxIterations, thermalConfig.relativeTolerance));
    heatTransferSolver->initialize(mesh, problemModel, materials);

    // Configure solid mechanics solver
    int mechOrder = getFieldOrder(FieldKind::Displacement);
    SolverConfig mechConfig = getSolverConfig(FieldKind::Displacement);
    solidMechanicsSolver->setOrder(mechOrder);
    solidMechanicsSolver->setSolver(createLinearSolver(
        mechConfig.type, mechConfig.maxIterations, mechConfig.relativeTolerance));
    solidMechanicsSolver->initialize(mesh, problemModel, materials);

    // Create coupling manager
    CouplingManager couplingManager;
    couplingManager.registerField(FieldKind::ElectricPotential, electrostaticsSolver.get());
    couplingManager.registerField(FieldKind::Temperature, heatTransferSolver.get());
    couplingManager.registerField(FieldKind::Displacement, solidMechanicsSolver.get());
    couplingManager.setCouplingConfig(problemModel.couplingConfig);

    // Register coupling types for dependency graph construction
    for (CouplingKind kind : problemModel.couplings) {
        couplingManager.registerCoupling(kind);
    }

    // Run coupled simulation
    couplingManager.run();

    // Extract results
    const mfem::GridFunction& potential = electrostaticsSolver->getField();
    const mfem::GridFunction& temperature = heatTransferSolver->getField();
    const mfem::GridFunction& displacement = solidMechanicsSolver->getField();

    recoverVertexSample(mesh, potential, temperature, displacement, result);
}

} // namespace mpfem
