#ifndef MPFEM_COUPLING_MANAGER_HPP
#define MPFEM_COUPLING_MANAGER_HPP

#include "physics_field_solver.hpp"
#include "physics_problem_model.hpp"

#include <map>
#include <memory>
#include <vector>

namespace mpfem {

/**
 * @brief Unified coupling manager for multi-physics simulations.
 * 
 * This class manages the coupling between different physics fields,
 * handling data transfer and iteration control. It supports:
 * - Joule heating coupling (electrostatics -> heat transfer)
 * - Thermal expansion coupling (heat transfer -> solid mechanics)
 * - Newton-Raphson and Picard iteration methods
 */
class CouplingManager {
public:
    CouplingManager();
    ~CouplingManager();

    /**
     * @brief Register a physics field solver.
     * @param kind The field kind this solver handles.
     * @param solver Pointer to the solver (ownership not transferred).
     */
    void registerField(FieldKind kind, PhysicsFieldSolver* solver);

    /**
     * @brief Set the coupling configuration.
     */
    void setCouplingConfig(const CouplingConfig& config);

    /**
     * @brief Setup coupling data transfer between fields.
     */
    void setupCoupling();

    /**
     * @brief Run the coupled iteration loop.
     */
    void run();

    /**
     * @brief Get a field solution by kind.
     */
    const mfem::GridFunction* getField(FieldKind kind) const;

    /**
     * @brief Get the number of iterations from the last run.
     */
    int getNumIterations() const;

    /**
     * @brief Get the convergence status from the last run.
     */
    bool isConverged() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace mpfem

#endif