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
    
    // Analyze coupling kinds to determine field dependencies
    bool hasJouleHeating = false;
    bool hasThermalExpansion = false;
    
    for (CouplingKind kind : couplingKinds_) {
        if (kind == CouplingKind::JouleHeating) {
            hasJouleHeating = true;
            // ElectricPotential <-> Temperature: bidirectional coupling
            tightlyCoupledFields_.insert(FieldKind::ElectricPotential);
            tightlyCoupledFields_.insert(FieldKind::Temperature);
        } else if (kind == CouplingKind::ThermalExpansion) {
            hasThermalExpansion = true;
            // Temperature -> Displacement: unidirectional dependency
            // Displacement is downstream, not part of tight coupling
        }
    }
    
    // Build downstream fields list (fields that depend on tightly coupled fields)
    if (hasThermalExpansion && solvers_.count(FieldKind::Displacement)) {
        downstreamFields_.push_back(FieldKind::Displacement);
    }
    
    Logger::log(LogLevel::Debug, "Dependency graph built: " 
               + std::to_string(tightlyCoupledFields_.size()) + " tightly coupled fields, "
               + std::to_string(downstreamFields_.size()) + " downstream fields");
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
        auto iter = solvers_.find(kind);
        if (iter != solvers_.end()) {
            const mfem::GridFunction& field = iter->second->getField();
            previousSolutions_.push_back(mfem::Vector(field));
        }
    }
}

bool CouplingManager::checkConvergence()
{
    if (previousSolutions_.empty()) {
        return false;
    }

    int index = 0;
    double maxRelDiff = 0.0;

    for (FieldKind kind : tightlyCoupledFields_) {
        auto iter = solvers_.find(kind);
        if (iter == solvers_.end()) continue;
        
        const mfem::GridFunction& current = iter->second->getField();
        const mfem::Vector& previous = previousSolutions_[index];
        
        double diff = 0.0;
        double norm = 0.0;
        const int size = current.Size();
        
        for (int i = 0; i < size; ++i) {
            const double d = current(i) - previous(i);
            diff += d * d;
            norm += current(i) * current(i);
        }

        const double relDiff = (norm > 0.0) ? std::sqrt(diff / norm) : 0.0;
        maxRelDiff = std::max(maxRelDiff, relDiff);
        ++index;
    }

    Logger::log(LogLevel::Debug, "Convergence check: max relative diff = " 
               + std::to_string(maxRelDiff));

    return maxRelDiff < config_.tolerance;
}

void CouplingManager::runTightlyCoupledIteration()
{
    // Step 1: Solve electrostatics (with temperature-dependent conductivity)
    auto electroIter = solvers_.find(FieldKind::ElectricPotential);
    if (electroIter != solvers_.end() && tightlyCoupledFields_.count(FieldKind::ElectricPotential)) {
        auto* electroSolver = dynamic_cast<ElectrostaticsSolver*>(electroIter->second);

        // Update temperature field for conductivity
        auto tempIter = solvers_.find(FieldKind::Temperature);
        if (tempIter != solvers_.end() && electroSolver) {
            electroSolver->setTemperatureField(&tempIter->second->getField());
        }

        solveField(electroIter->second, "Electrostatics solve");
    }

    // Step 2: Solve heat transfer (with Joule heating source)
    auto heatIter = solvers_.find(FieldKind::Temperature);
    if (heatIter != solvers_.end() && tightlyCoupledFields_.count(FieldKind::Temperature)) {
        auto* heatSolver = dynamic_cast<HeatTransferSolver*>(heatIter->second);
        Check(heatSolver != nullptr, "Failed to cast Temperature solver to HeatTransferSolver");

        // Update potential field and conductivity for Joule heating
        if (electroIter != solvers_.end()) {
            auto* electroSolver = dynamic_cast<ElectrostaticsSolver*>(electroIter->second);
            Check(electroSolver != nullptr, "Failed to cast ElectricPotential solver to ElectrostaticsSolver");
            
            mfem::Coefficient* condCoef = electroSolver->getConductivityCoefficient();
            Check(condCoef != nullptr, "Conductivity coefficient is null");
            
            heatSolver->setPotentialField(&electroIter->second->getField());
            heatSolver->setConductivityCoefficient(condCoef);
        }

        solveField(heatIter->second, "Heat transfer solve");
    }
}

void CouplingManager::solveDownstreamFields()
{
    // Solve solid mechanics (with thermal expansion load) - only once after convergence
    auto heatIter = solvers_.find(FieldKind::Temperature);
    
    for (FieldKind kind : downstreamFields_) {
        if (kind == FieldKind::Displacement) {
            auto mechIter = solvers_.find(FieldKind::Displacement);
            if (mechIter != solvers_.end()) {
                auto* mechSolver = dynamic_cast<SolidMechanicsSolver*>(mechIter->second);

                // Update temperature field for thermal expansion
                if (heatIter != solvers_.end() && mechSolver) {
                    mechSolver->setTemperatureField(&heatIter->second->getField());
                }

                solveField(mechIter->second, "Solid mechanics solve");
            }
        }
    }
}

void CouplingManager::run()
{
    converged_ = false;
    numIterations_ = 0;

    // Build dependency graph based on registered couplings
    buildDependencyGraph();

    ScopedTimer timer("Coupled solve");

    // Phase 1: Iterate tightly coupled fields until convergence
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
    auto iter = solvers_.find(kind);
    if (iter != solvers_.end()) {
        return &iter->second->getField();
    }
    return nullptr;
}

} // namespace mpfem
