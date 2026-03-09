#include "coupling_manager.hpp"
#include "electrostatics_solver.hpp"
#include "heat_transfer_solver.hpp"
#include "solid_mechanics_solver.hpp"
#include "logger.hpp"

#include <cmath>

namespace mpfem {

class CouplingManager::Impl {
public:
    std::map<FieldKind, PhysicsFieldSolver*> solvers_;
    CouplingConfig config_;
    int numIterations_ = 0;
    bool converged_ = false;

    void runNewtonRaphson();
    void runPicard();
    void solveField(PhysicsFieldSolver* solver, const std::string& name);
    void savePreviousSolutions();
    bool checkConvergence();

    std::vector<mfem::Vector> previousSolutions_;
};

void CouplingManager::Impl::solveField(PhysicsFieldSolver* solver, const std::string& name)
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

void CouplingManager::Impl::savePreviousSolutions()
{
    previousSolutions_.clear();
    for (const auto& [kind, solver] : solvers_) {
        const mfem::GridFunction& field = solver->getField();
        previousSolutions_.push_back(mfem::Vector(field));
    }
}

bool CouplingManager::Impl::checkConvergence()
{
    if (previousSolutions_.empty()) {
        return false;
    }

    int index = 0;
    double maxRelDiff = 0.0;

    for (const auto& [kind, solver] : solvers_) {
        const mfem::GridFunction& current = solver->getField();
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

void CouplingManager::Impl::runNewtonRaphson()
{
    for (int iter = 0; iter < config_.maxIterations; ++iter) {
        numIterations_ = iter + 1;
        savePreviousSolutions();

        Logger::log(LogLevel::Info, "--- Coupling iteration " + std::to_string(iter + 1) + " ---");

        // Step 1: Solve electrostatics (with temperature-dependent conductivity)
        auto electroIter = solvers_.find(FieldKind::ElectricPotential);
        if (electroIter != solvers_.end()) {
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
        if (heatIter != solvers_.end()) {
            auto* heatSolver = dynamic_cast<HeatTransferSolver*>(heatIter->second);

            // Update potential field and conductivity for Joule heating
            if (electroIter != solvers_.end() && heatSolver) {
                heatSolver->setPotentialField(&electroIter->second->getField());
                heatSolver->setConductivityCoefficient(
                    dynamic_cast<ElectrostaticsSolver*>(electroIter->second)
                        ->getConductivityCoefficient());
            }

            solveField(heatIter->second, "Heat transfer solve");
        }

        // Step 3: Solve solid mechanics (with thermal expansion load)
        auto mechIter = solvers_.find(FieldKind::Displacement);
        if (mechIter != solvers_.end()) {
            auto* mechSolver = dynamic_cast<SolidMechanicsSolver*>(mechIter->second);

            // Update temperature field for thermal expansion
            if (heatIter != solvers_.end() && mechSolver) {
                mechSolver->setTemperatureField(&heatIter->second->getField());
            }

            solveField(mechIter->second, "Solid mechanics solve");
        }

        // Check convergence
        if (checkConvergence()) {
            converged_ = true;
            Logger::log(LogLevel::Info, "Coupling converged after " + std::to_string(numIterations_) + " iterations");
            return;
        }
    }

    Logger::log(LogLevel::Error, "Coupling did not converge after " 
               + std::to_string(config_.maxIterations) + " iterations");
    converged_ = false;
}

void CouplingManager::Impl::runPicard()
{
    runNewtonRaphson();
}

CouplingManager::CouplingManager()
    : impl_(std::make_unique<Impl>())
{
}

CouplingManager::~CouplingManager() = default;

void CouplingManager::registerField(FieldKind kind, PhysicsFieldSolver* solver)
{
    impl_->solvers_[kind] = solver;
}

void CouplingManager::setCouplingConfig(const CouplingConfig& config)
{
    impl_->config_ = config;
}

void CouplingManager::setupCoupling()
{
    // Setup coupling connections between fields
}

void CouplingManager::run()
{
    impl_->converged_ = false;
    impl_->numIterations_ = 0;

    ScopedTimer timer("Coupled solve");

    CouplingMethod method = parseCouplingMethod(impl_->config_.method);
    if (method == CouplingMethod::Picard) {
        impl_->runPicard();
    } else {
        impl_->runNewtonRaphson();
    }
}

const mfem::GridFunction* CouplingManager::getField(FieldKind kind) const
{
    auto iter = impl_->solvers_.find(kind);
    if (iter != impl_->solvers_.end()) {
        return &iter->second->getField();
    }
    return nullptr;
}

int CouplingManager::getNumIterations() const
{
    return impl_->numIterations_;
}

bool CouplingManager::isConverged() const
{
    return impl_->converged_;
}

} // namespace mpfem