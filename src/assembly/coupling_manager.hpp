#ifndef MPFEM_COUPLING_MANAGER_HPP
#define MPFEM_COUPLING_MANAGER_HPP

#include "physics_field_solver.hpp"
#include "physics_problem_model.hpp"

#include <map>
#include <memory>
#include <vector>
#include <set>

namespace mpfem {

/**
 * @brief Unified coupling manager for multi-physics simulations.
 * 
 * Optimized coupling strategy:
 * 1. Build dependency graph from CouplingKind definitions
 * 2. Identify tightly coupled subgraphs (bidirectional dependencies)
 * 3. Iterate only tightly coupled fields until convergence
 * 4. Solve downstream fields (unidirectional dependencies) once after convergence
 */
class CouplingManager {
public:
    CouplingManager() = default;
    ~CouplingManager() = default;

    /**
     * @brief Register a physics field solver.
     */
    void registerField(FieldKind kind, PhysicsFieldSolver* solver);

    /**
     * @brief Set the coupling configuration.
     */
    void setCouplingConfig(const CouplingConfig& config);

    /**
     * @brief Register a coupling type between fields.
     */
    void registerCoupling(CouplingKind kind);

    /**
     * @brief Run the coupled iteration loop with optimized dependency resolution.
     */
    void run();

    /**
     * @brief Get a field solution by kind.
     */
    const mfem::GridFunction* getField(FieldKind kind) const;

    /**
     * @brief Get the number of iterations from the last run.
     */
    int getNumIterations() const { return numIterations_; }

    /**
     * @brief Get the convergence status from the last run.
     */
    bool isConverged() const { return converged_; }

private:
    // Field management
    std::map<FieldKind, PhysicsFieldSolver*> solvers_;
    CouplingConfig config_;
    
    // Iteration state
    int numIterations_ = 0;
    bool converged_ = false;
    std::vector<mfem::Vector> previousSolutions_;
    
    // Dependency graph
    std::set<FieldKind> tightlyCoupledFields_;  // Fields that need iteration
    std::vector<FieldKind> downstreamFields_;   // Fields solved after convergence
    std::vector<CouplingKind> couplingKinds_;
    
    // Build dependency graph from coupling kinds
    void buildDependencyGraph();
    
    // Solve a single field
    void solveField(PhysicsFieldSolver* solver, const std::string& name);
    
    // Save solutions for convergence check
    void savePreviousSolutions();
    
    // Check convergence for tightly coupled fields only
    bool checkConvergence();
    
    // Run iteration for tightly coupled fields
    void runTightlyCoupledIteration();
    
    // Solve downstream fields after convergence
    void solveDownstreamFields();
};

} // namespace mpfem

#endif
