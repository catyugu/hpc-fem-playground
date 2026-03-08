#include "coupling_manager.hpp"
#include "electrostatics_solver.hpp"
#include "heat_transfer_solver.hpp"
#include "solid_mechanics_solver.hpp"

#include <cmath>
#include <algorithm>

namespace mpfem {

class CouplingManager::Impl {
public:
    std::map<FieldKind, PhysicsFieldSolver*> solvers_;
    CouplingConfig config_;
    int numIterations_ = 0;
    bool converged_ = false;

    // Previous iteration solutions for convergence check
    std::unique_ptr<mfem::GridFunction> prevTemperature_;
    std::unique_ptr<mfem::GridFunction> prevDisplacement_;

    double computeFieldDelta(const mfem::GridFunction& previous,
                            const mfem::GridFunction& current);
    bool runNewtonRaphson(std::string& errorMessage);
    bool runPicard(std::string& errorMessage);
    void savePreviousSolutions();
    bool checkConvergence();
};

double CouplingManager::Impl::computeFieldDelta(const mfem::GridFunction& previous,
                                                 const mfem::GridFunction& current)
{
    if (previous.Size() != current.Size()) {
        return 1.0;
    }

    double maxDelta = 0.0;
    for (int i = 0; i < previous.Size(); ++i) {
        const double delta = std::abs(previous(i) - current(i));
        if (delta > maxDelta) {
            maxDelta = delta;
        }
    }

    return maxDelta;
}

void CouplingManager::Impl::savePreviousSolutions()
{
    // Save temperature
    auto tempIter = solvers_.find(FieldKind::Temperature);
    if (tempIter != solvers_.end()) {
        const auto& temp = tempIter->second->getField();
        prevTemperature_ = std::make_unique<mfem::GridFunction>(temp);
    }

    // Save displacement
    auto dispIter = solvers_.find(FieldKind::Displacement);
    if (dispIter != solvers_.end()) {
        const auto& disp = dispIter->second->getField();
        prevDisplacement_ = std::make_unique<mfem::GridFunction>(disp);
    }
}

bool CouplingManager::Impl::checkConvergence()
{
    // Check temperature convergence
    auto tempIter = solvers_.find(FieldKind::Temperature);
    if (tempIter != solvers_.end() && prevTemperature_) {
        const auto& currentTemp = tempIter->second->getField();
        const double tempDelta = computeFieldDelta(*prevTemperature_, currentTemp);
        if (tempDelta >= config_.tolerance) {
            return false;
        }
    }

    // Check displacement convergence
    auto dispIter = solvers_.find(FieldKind::Displacement);
    if (dispIter != solvers_.end() && prevDisplacement_) {
        const auto& currentDisp = dispIter->second->getField();
        const double dispDelta = computeFieldDelta(*prevDisplacement_, currentDisp);
        if (dispDelta >= config_.tolerance) {
            return false;
        }
    }

    return true;
}

bool CouplingManager::Impl::runNewtonRaphson(std::string& errorMessage)
{
    // For segregated coupling with Newton-Raphson, we use a simplified approach:
    // Each field is solved sequentially with updated coupling terms.
    // This is equivalent to a block Gauss-Seidel iteration with under-relaxation.

    for (int iter = 0; iter < config_.maxIterations; ++iter) {
        numIterations_ = iter + 1;
        savePreviousSolutions();

        // Step 1: Solve electrostatics (with temperature-dependent conductivity)
        auto electroIter = solvers_.find(FieldKind::ElectricPotential);
        if (electroIter != solvers_.end()) {
            auto* electroSolver = dynamic_cast<ElectrostaticsSolver*>(electroIter->second);

            // Update temperature field for conductivity
            auto tempIter = solvers_.find(FieldKind::Temperature);
            if (tempIter != solvers_.end() && electroSolver) {
                electroSolver->setTemperatureField(&tempIter->second->getField());
            }

            // Apply boundary conditions before each solve
            electroIter->second->applyBoundaryConditions();

            if (!electroIter->second->assemble(errorMessage)) {
                return false;
            }
            if (!electroIter->second->solve(errorMessage)) {
                return false;
            }
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

            // Apply boundary conditions before each solve
            heatIter->second->applyBoundaryConditions();

            if (!heatIter->second->assemble(errorMessage)) {
                return false;
            }
            if (!heatIter->second->solve(errorMessage)) {
                return false;
            }
        }

        // Step 3: Solve solid mechanics (with thermal expansion load)
        auto mechIter = solvers_.find(FieldKind::Displacement);
        if (mechIter != solvers_.end()) {
            auto* mechSolver = dynamic_cast<SolidMechanicsSolver*>(mechIter->second);

            // Update temperature field for thermal expansion
            if (heatIter != solvers_.end() && mechSolver) {
                mechSolver->setTemperatureField(&heatIter->second->getField());
            }

            // Apply boundary conditions before each solve
            mechIter->second->applyBoundaryConditions();

            if (!mechIter->second->assemble(errorMessage)) {
                return false;
            }
            if (!mechIter->second->solve(errorMessage)) {
                return false;
            }
        }

        // Check convergence
        if (checkConvergence()) {
            converged_ = true;
            return true;
        }
    }

    errorMessage = "Newton-Raphson coupling did not converge after "
                 + std::to_string(config_.maxIterations) + " iterations";
    converged_ = false;
    return false;
}

bool CouplingManager::Impl::runPicard(std::string& errorMessage)
{
    // Picard iteration is similar to Newton-Raphson for segregated coupling
    // The difference is mainly in how we handle the linearization
    // For this implementation, we use the same approach as Newton-Raphson
    return runNewtonRaphson(errorMessage);
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
    // Initial setup - nothing specific needed here for segregated coupling
}

bool CouplingManager::run(std::string& errorMessage)
{
    errorMessage.clear();
    impl_->converged_ = false;
    impl_->numIterations_ = 0;

    // Initial boundary condition application
    for (auto& pair : impl_->solvers_) {
        pair.second->applyBoundaryConditions();
    }

    // Choose iteration method
    if (impl_->config_.method == "newton_raphson") {
        return impl_->runNewtonRaphson(errorMessage);
    } else {
        return impl_->runPicard(errorMessage);
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
