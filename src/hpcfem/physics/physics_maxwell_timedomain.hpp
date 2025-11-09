/**
 * @file physics_maxwell_timedomain.hpp
 * @brief Time-domain Maxwell equation solver with mixed formulation
 * 
 * Solves the coupled first-order Maxwell equations:
 *   ε ∂E/∂t = ∇×(μ⁻¹B) - σE - J
 *   ∂B/∂t = -∇×E
 * 
 * Uses H(curl) discretization for E and H(div) for B with energy-conserving
 * symplectic time integration. Supports lossy materials, sources, and
 * sophisticated boundary conditions.
 */

#ifndef HPCFEM_PHYSICS_MAXWELL_TIMEDOMAIN_HPP
#define HPCFEM_PHYSICS_MAXWELL_TIMEDOMAIN_HPP

#include "mfem.hpp"
#include "maxwell_materials.hpp"
#include "maxwell_sources.hpp"
#include "maxwell_boundary.hpp"

namespace hpcfem
{

/**
 * @brief Time-domain Maxwell solver using mixed formulation
 * 
 * This class implements a time-dependent operator for the coupled
 * first-order Maxwell equations. It inherits from MFEM's
 * TimeDependentOperator to enable use with various time integrators.
 * 
 * The operator implements:
 *   dy/dt = f(t, y)
 * where y = [E; B] contains both E-field and B-field DOFs.
 */
class PhysicsMaxwellTimeDomain : public mfem::TimeDependentOperator
{
public:
    /**
     * @brief Constructor for parallel execution
     * @param pmesh Parallel mesh
     * @param order Polynomial order for finite elements
     * @param epsilonFunc Electric permittivity function ε(x)
     * @param muInvFunc Inverse magnetic permeability μ⁻¹(x)
     * @param sigmaFunc Electric conductivity σ(x), can be nullptr
     * @param jFunc Volumetric current density J(x,t), can be nullptr
     * @param abcMarkers Boundary attributes for absorbing BC
     * @param dbcMarkers Boundary attributes for Dirichlet BC
     * @param dEdtFunc Time derivative of E on Dirichlet boundaries
     */
#ifdef MFEM_USE_MPI
    PhysicsMaxwellTimeDomain(mfem::ParMesh* pmesh,
                            int order,
                            double (*epsilonFunc)(const mfem::Vector&),
                            double (*muInvFunc)(const mfem::Vector&),
                            double (*sigmaFunc)(const mfem::Vector&),
                            void (*jFunc)(const mfem::Vector&, double, mfem::Vector&),
                            const mfem::Array<int>& abcMarkers,
                            const mfem::Array<int>& dbcMarkers,
                            void (*dEdtFunc)(const mfem::Vector&, double, mfem::Vector&));
#else
    PhysicsMaxwellTimeDomain(mfem::Mesh* mesh,
                            int order,
                            double (*epsilonFunc)(const mfem::Vector&),
                            double (*muInvFunc)(const mfem::Vector&),
                            double (*sigmaFunc)(const mfem::Vector&),
                            void (*jFunc)(const mfem::Vector&, double, mfem::Vector&),
                            const mfem::Array<int>& abcMarkers,
                            const mfem::Array<int>& dbcMarkers,
                            void (*dEdtFunc)(const mfem::Vector&, double, mfem::Vector&));
#endif

    /**
     * @brief Destructor
     */
    ~PhysicsMaxwellTimeDomain();

    /**
     * @brief Compute dy/dt = f(t,y) for explicit time integration
     * @param B Magnetic flux vector (input)
     * @param dEdt Time derivative of E-field (output)
     */
    virtual void Mult(const mfem::Vector& B, mfem::Vector& dEdt) const;

    /**
     * @brief Implicit solve for lossy materials
     * @param dt Time step size
     * @param B Magnetic flux at time t
     * @param dEdt Time derivative of E at time t+dt (output)
     */
    void implicitSolve(double dt, const mfem::Vector& B, mfem::Vector& dEdt);

    /**
     * @brief Get maximum stable time step (CFL condition)
     * @return Maximum time step for stability
     */
    double getMaximumTimeStep() const;

    /**
     * @brief Compute total electromagnetic energy
     * @return Energy = 0.5*∫(ε|E|² + μ|B|²) dV
     */
    double getEnergy() const;

    /**
     * @brief Get E-field grid function
     */
#ifdef MFEM_USE_MPI
    mfem::ParGridFunction* getEField() { return eField_; }
#else
    mfem::GridFunction* getEField() { return eField_; }
#endif

    /**
     * @brief Get B-field grid function
     */
#ifdef MFEM_USE_MPI
    mfem::ParGridFunction* getBField() { return bField_; }
#else
    mfem::GridFunction* getBField() { return bField_; }
#endif

    /**
     * @brief Set initial E-field from coefficient
     */
    void setInitialEField(mfem::VectorCoefficient& eCoef);

    /**
     * @brief Set initial B-field from coefficient
     */
    void setInitialBField(mfem::VectorCoefficient& bCoef);

    /**
     * @brief Get E-field DOF vector
     */
#ifdef MFEM_USE_MPI
    mfem::HypreParVector* getEVector() { return eVec_; }
#else
    mfem::Vector* getEVector() { return eField_; }
#endif

    /**
     * @brief Get B-field DOF vector
     */
#ifdef MFEM_USE_MPI
    mfem::HypreParVector* getBVector() { return bVec_; }
#else
    mfem::Vector* getBVector() { return bField_; }
#endif

    /**
     * @brief Get negative curl operator for dB/dt = -∇×E
     */
#ifdef MFEM_USE_MPI
    mfem::HypreParMatrix* getNegCurlOperator() { return negCurlMatrix_; }
#else
    mfem::SparseMatrix* getNegCurlOperator() { return negCurlMatrix_; }
#endif

    /**
     * @brief Sync grid functions from DOF vectors
     */
    void syncGridFunctions();

private:
    /**
     * @brief Assemble system matrices
     */
    void assembleMatrices();

    /**
     * @brief Setup linear solver for implicit integration
     */
    void setupImplicitSolver(double dt);

    int order_;
    bool isLossy_;

#ifdef MFEM_USE_MPI
    mfem::ParMesh* pmesh_;
    mfem::ParFiniteElementSpace* hCurlFESpace_;
    mfem::ParFiniteElementSpace* hDivFESpace_;
    mfem::ParGridFunction* eField_;
    mfem::ParGridFunction* bField_;
    mfem::ParGridFunction* dEdtField_;
    mfem::ParBilinearForm* hDivMassMuInv_;
    mfem::ParBilinearForm* hCurlLosses_;
    mfem::ParMixedBilinearForm* weakCurlMuInv_;
    mfem::ParDiscreteLinearOperator* curlOp_;
    mfem::HypreParMatrix* massMuInvMatrix_;
    mfem::HypreParMatrix* lossMatrix_;
    mfem::HypreParMatrix* negCurlMatrix_;
    mfem::HypreParMatrix* weakCurlMatrix_;
    mfem::HypreParVector* eVec_;
    mfem::HypreParVector* bVec_;
    mfem::HypreParVector* rhsVec_;
#else
    mfem::Mesh* mesh_;
    mfem::FiniteElementSpace* hCurlFESpace_;
    mfem::FiniteElementSpace* hDivFESpace_;
    mfem::GridFunction* eField_;
    mfem::GridFunction* bField_;
    mfem::GridFunction* dEdtField_;
    mfem::BilinearForm* hDivMassMuInv_;
    mfem::BilinearForm* hCurlLosses_;
    mfem::MixedBilinearForm* weakCurlMuInv_;
    mfem::DiscreteLinearOperator* curlOp_;
    mfem::SparseMatrix* massMuInvMatrix_;
    mfem::SparseMatrix* lossMatrix_;
    mfem::SparseMatrix* negCurlMatrix_;
    mfem::SparseMatrix* weakCurlMatrix_;
#endif

    mfem::FiniteElementCollection* hCurlFEC_;
    mfem::FiniteElementCollection* hDivFEC_;

    MaxwellMaterialProperties* materials_;
    MaxwellSourceTerms* sources_;
    MaxwellBoundaryConditions* boundaries_;

    // Linear solver for implicit integration
    mfem::CGSolver* implicitSolver_;
    mfem::HypreSmoother* preconditioner_;
    
    double maxTimeStep_;
    
    // MPI info
    int myRank_;
    int numProcs_;
};

} // namespace hpcfem

#endif // HPCFEM_PHYSICS_MAXWELL_TIMEDOMAIN_HPP
