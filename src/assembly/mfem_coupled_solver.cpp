#include "mfem_coupled_solver.hpp"

#include "coupling_manager.hpp"
#include "electrostatics_solver.hpp"
#include "heat_transfer_solver.hpp"
#include "linear_solver_strategy.hpp"
#include "solid_mechanics_solver.hpp"
#include "logger.hpp"
#include "mpfem_types.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <map>
#include <vector>

namespace mpfem {

namespace {

const int RESULT_COMPONENT_COUNT = 6;
const double COORDINATE_KEY_SCALE = 1.0e12;

struct CoordinateKey {
    std::int64_t x = 0;
    std::int64_t y = 0;
    std::int64_t z = 0;
};

bool operator<(const CoordinateKey& left, const CoordinateKey& right)
{
    if (left.x != right.x) {
        return left.x < right.x;
    }
    if (left.y != right.y) {
        return left.y < right.y;
    }
    return left.z < right.z;
}

CoordinateKey makeCoordinateKey(double x, double y, double z)
{
    CoordinateKey key;
    key.x = static_cast<std::int64_t>(std::llround(x * COORDINATE_KEY_SCALE));
    key.y = static_cast<std::int64_t>(std::llround(y * COORDINATE_KEY_SCALE));
    key.z = static_cast<std::int64_t>(std::llround(z * COORDINATE_KEY_SCALE));
    return key;
}

void recoverVertexSampleSerial(const FemMesh& mesh,
                               const FemGridFunction& potential,
                               const FemGridFunction& temperature,
                               const FemGridFunction& displacement,
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

void recoverVertexSample(const FemMesh& mesh,
                         const FemGridFunction& potential,
                         const FemGridFunction& temperature,
                         const FemGridFunction& displacement,
                         CoupledFieldResult& result)
{
#ifdef MFEM_USE_MPI
    const int localVertexCount = mesh.GetNV();

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

    Check(nodalPotential.Size() >= localVertexCount, "Nodal potential size mismatch");
    Check(nodalTemperature.Size() >= localVertexCount, "Nodal temperature size mismatch");
    Check(nodalUx.Size() >= localVertexCount, "Nodal Ux size mismatch");
    Check(nodalUy.Size() >= localVertexCount, "Nodal Uy size mismatch");
    Check(nodalUz.Size() >= localVertexCount, "Nodal Uz size mismatch");

    std::vector<double> localPacked(localVertexCount * RESULT_COMPONENT_COUNT, 0.0);
    for (int i = 0; i < localVertexCount; ++i) {
        const double* vertex = mesh.GetVertex(i);
        const int offset = i * RESULT_COMPONENT_COUNT;
        localPacked[offset + 0] = vertex[0];
        localPacked[offset + 1] = vertex[1];
        localPacked[offset + 2] = vertex[2];
        localPacked[offset + 3] = nodalPotential(i);
        localPacked[offset + 4] = nodalTemperature(i);

        const double ux = nodalUx(i);
        const double uy = nodalUy(i);
        const double uz = nodalUz(i);
        localPacked[offset + 5] = std::sqrt(ux * ux + uy * uy + uz * uz);
    }

    const int rank = mesh.GetMyRank();
    const int size = mesh.GetNRanks();
    std::vector<int> recvCounts;
    if (rank == 0) {
        recvCounts.resize(size, 0);
    }
    int* recvCountsPtr = rank == 0 ? recvCounts.data() : nullptr;

    const int localCount = static_cast<int>(localPacked.size());
    MPI_Gather(&localCount, 1, MPI_INT,
               recvCountsPtr, 1, MPI_INT, 0, mesh.GetComm());

    std::vector<int> displacements;
    std::vector<double> gathered;
    if (rank == 0) {
        displacements.resize(size, 0);
        for (int i = 1; i < size; ++i) {
            displacements[i] = displacements[i - 1] + recvCounts[i - 1];
        }
        const int totalCount = displacements[size - 1] + recvCounts[size - 1];
        gathered.resize(totalCount, 0.0);
    }
    int* displacementsPtr = rank == 0 ? displacements.data() : nullptr;
    double* gatheredPtr = rank == 0 ? gathered.data() : nullptr;

    MPI_Gatherv(localPacked.data(), localCount, MPI_DOUBLE,
                gatheredPtr, recvCountsPtr, displacementsPtr, MPI_DOUBLE,
                0, mesh.GetComm());

    if (rank != 0) {
        result = CoupledFieldResult();
        return;
    }

    std::map<CoordinateKey, std::array<double, RESULT_COMPONENT_COUNT>> uniqueSamples;
    for (std::size_t i = 0; i + RESULT_COMPONENT_COUNT <= gathered.size(); i += RESULT_COMPONENT_COUNT) {
        const double x = gathered[i + 0];
        const double y = gathered[i + 1];
        const double z = gathered[i + 2];
        const CoordinateKey key = makeCoordinateKey(x, y, z);

        std::array<double, RESULT_COMPONENT_COUNT> values;
        values[0] = x;
        values[1] = y;
        values[2] = z;
        values[3] = gathered[i + 3];
        values[4] = gathered[i + 4];
        values[5] = gathered[i + 5];
        uniqueSamples[key] = values;
    }

    const std::size_t pointCount = uniqueSamples.size();
    result.coordinates.resize(pointCount);
    result.electricPotential.resize(pointCount, 0.0);
    result.temperature.resize(pointCount, 0.0);
    result.displacement.resize(pointCount, 0.0);

    std::size_t index = 0;
    for (const auto& entry : uniqueSamples) {
        const std::array<double, RESULT_COMPONENT_COUNT>& values = entry.second;
        result.coordinates[index].x = values[0];
        result.coordinates[index].y = values[1];
        result.coordinates[index].z = values[2];
        result.electricPotential[index] = values[3];
        result.temperature[index] = values[4];
        result.displacement[index] = values[5];
        ++index;
    }
#else
    recoverVertexSampleSerial(mesh, potential, temperature, displacement, result);
#endif
}

} // namespace

void MfemCoupledSolver::solve(FemMesh& mesh,
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

    // Setup coupling
    couplingManager.setupCoupling();

    // Run coupled simulation
    couplingManager.run();

    // Extract results
    const FemGridFunction& potential = electrostaticsSolver->getField();
    const FemGridFunction& temperature = heatTransferSolver->getField();
    const FemGridFunction& displacement = solidMechanicsSolver->getField();

    recoverVertexSample(mesh, potential, temperature, displacement, result);
}

} // namespace mpfem
