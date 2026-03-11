#include "coupling_manager.hpp"
#include "electrostatics_solver.hpp"
#include "heat_transfer_solver.hpp"
#include "solid_mechanics_solver.hpp"
#include "logger.hpp"

#include <cmath>

namespace mpfem {

void CouplingManager::registerField(FieldKind kind, PhysicsFieldSolver* solver)
{
    solvers_[kind] = solver;
}

void CouplingManager::setCouplingConfig(const CouplingConfig& config)
{
    config_ = config;
}

void CouplingManager::registerCoupling(CouplingKind kind)
{
    couplingKinds_.push_back(kind);
}

void CouplingManager::buildDependencyGraph()
{
    tightlyCoupledFields_.clear();
    downstreamFields_.clear();

    bool hasJouleHeating = false;
    bool hasThermalExpansion = false;

    for (CouplingKind kind : couplingKinds_) {
        if (kind == CouplingKind::JouleHeating) {
            hasJouleHeating = true;
            tightlyCoupledFields_.insert(FieldKind::ElectricPotential);
            tightlyCoupledFields_.insert(FieldKind::Temperature);
        } else if (kind == CouplingKind::ThermalExpansion) {
            hasThermalExpansion = true;
        }
    }

    if (hasThermalExpansion && solvers_.count(FieldKind::Displacement)) {
        downstreamFields_.push_back(FieldKind::Displacement);
    }

    Logger::log(LogLevel::Debug, "Dependency graph: "
               + std::to_string(tightlyCoupledFields_.size()) + " tightly coupled, "
               + std::to_string(downstreamFields_.size()) + " downstream");
}

void CouplingManager::solveField(PhysicsFieldSolver* solver, const std::string& name)
{
    ScopedTimer totalTimer(name);

    ScopedTimer bcTimer(name + " BC");
    solver->applyBoundaryConditions();
    bcTimer.stop();

    ScopedTimer assembleTimer(name + " assemble");
    solver->assemble();
    assembleTimer.stop();

    ScopedTimer solveTimer(name + " solve");
    solver->solve();
    solveTimer.stop();

    Logger::log(LogLevel::Info, name + " breakdown: BC=" + bcTimer.getElapsedStr()
        + ", assemble=" + assembleTimer.getElapsedStr()
        + ", solve=" + solveTimer.getElapsedStr());
}

void CouplingManager::savePreviousSolutions()
{
    previousSolutions_.clear();
    for (FieldKind kind : tightlyCoupledFields_) {
        auto it = solvers_.find(kind);
        if (it != solvers_.end()) {
            previousSolutions_.push_back(mfem::Vector(it->second->getField()));
        }
    }
}

bool CouplingManager::checkConvergence()
{
    if (previousSolutions_.empty()) return false;

    int idx = 0;
    double maxRelDiff = 0.0;

    for (FieldKind kind : tightlyCoupledFields_) {
        auto it = solvers_.find(kind);
        if (it == solvers_.end()) continue;

        const mfem::GridFunction& current = it->second->getField();
        const mfem::Vector& prev = previousSolutions_[idx];

        double diff = 0.0, norm = 0.0;
        const int n = current.Size();
        for (int i = 0; i < n; ++i) {
            const double d = current(i) - prev(i);
            diff += d * d;
            norm += current(i) * current(i);
        }

        const double relDiff = (norm > 0.0) ? std::sqrt(diff / norm) : 0.0;
        maxRelDiff = std::max(maxRelDiff, relDiff);
        ++idx;
    }

    Logger::log(LogLevel::Debug, "Convergence: max rel diff = " + std::to_string(maxRelDiff));
    return maxRelDiff < config_.tolerance;
}

void CouplingManager::runTightlyCoupledIteration()
{
    // Cache solver pointers (only do lookup once)
    auto electroIt = solvers_.find(FieldKind::ElectricPotential);
    auto heatIt = solvers_.find(FieldKind::Temperature);

    ElectrostaticsSolver* electroSolver = nullptr;
    HeatTransferSolver* heatSolver = nullptr;

    if (electroIt != solvers_.end() && tightlyCoupledFields_.count(FieldKind::ElectricPotential)) {
        electroSolver = dynamic_cast<ElectrostaticsSolver*>(electroIt->second);
    }
    if (heatIt != solvers_.end() && tightlyCoupledFields_.count(FieldKind::Temperature)) {
        heatSolver = dynamic_cast<HeatTransferSolver*>(heatIt->second);
        Check(heatSolver != nullptr, "Temperature solver is not a HeatTransferSolver");
    }

    // Step 1: Solve electrostatics (with temperature-dependent conductivity)
    if (electroSolver && heatIt != solvers_.end()) {
        electroSolver->setTemperatureField(&heatIt->second->getField());
        solveField(electroIt->second, "Electrostatics solve");
    }

    // Step 2: Solve heat transfer (with Joule heating source)
    if (heatSolver && electroIt != solvers_.end()) {
        Check(electroSolver != nullptr, "ElectricPotential solver is not an ElectrostaticsSolver");

        mfem::Coefficient* condCoef = electroSolver->getConductivityCoefficient();
        Check(condCoef != nullptr, "Conductivity coefficient is null");

        heatSolver->setPotentialField(&electroIt->second->getField());
        heatSolver->setConductivityCoefficient(condCoef);
        solveField(heatIt->second, "Heat transfer solve");
    }
}

void CouplingManager::solveDownstreamFields()
{
    auto heatIt = solvers_.find(FieldKind::Temperature);

    for (FieldKind kind : downstreamFields_) {
        if (kind == FieldKind::Displacement) {
            auto mechIt = solvers_.find(FieldKind::Displacement);
            if (mechIt != solvers_.end()) {
                auto* mechSolver = dynamic_cast<SolidMechanicsSolver*>(mechIt->second);
                if (mechSolver && heatIt != solvers_.end()) {
                    mechSolver->setTemperatureField(&heatIt->second->getField());
                }
                solveField(mechIt->second, "Solid mechanics solve");
            }
        }
    }
}

void CouplingManager::run()
{
    converged_ = false;
    numIterations_ = 0;

    buildDependencyGraph();

    ScopedTimer timer("Coupled solve");

    // Phase 1: Iterate tightly coupled fields
    Logger::log(LogLevel::Info, "=== Phase 1: Tightly coupled iteration ===");

    for (int iter = 0; iter < config_.maxIterations; ++iter) {
        numIterations_ = iter + 1;
        savePreviousSolutions();

        Logger::log(LogLevel::Info, "--- Coupling iteration " + std::to_string(iter + 1) + " ---");

        runTightlyCoupledIteration();

        if (checkConvergence()) {
            converged_ = true;
            Logger::log(LogLevel::Info, "Tightly coupled fields converged after "
                       + std::to_string(numIterations_) + " iterations");
            break;
        }
    }

    if (!converged_) {
        Logger::log(LogLevel::Error, "Coupling did not converge after "
                   + std::to_string(config_.maxIterations) + " iterations");
    }

    // Phase 2: Solve downstream fields once
    if (!downstreamFields_.empty()) {
        Logger::log(LogLevel::Info, "=== Phase 2: Solving downstream fields ===");
        solveDownstreamFields();
    }
}

const mfem::GridFunction* CouplingManager::getField(FieldKind kind) const
{
    auto it = solvers_.find(kind);
    return (it != solvers_.end()) ? &it->second->getField() : nullptr;
}

} // namespace mpfem