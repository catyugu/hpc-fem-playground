/**
 * @file maxwell_boundary.hpp
 * @brief Boundary condition management for Maxwell time-domain solver
 * 
 * Handles various boundary conditions:
 * - Natural (zero tangential current)
 * - Dirichlet (prescribed E-field time derivative)
 * - Sommerfeld absorbing boundary condition (ABC)
 * 
 * TODO: PML (Perfectly Matched Layer)
 * TODO: Port boundary conditions for waveguides
 */

#ifndef HPCFEM_MAXWELL_BOUNDARY_HPP
#define HPCFEM_MAXWELL_BOUNDARY_HPP

#include "mfem.hpp"

namespace hpcfem
{

/**
 * @brief Boundary condition types for Maxwell solver
 */
enum class MaxwellBoundaryType
{
    Natural,    // Zero tangential current (default, do nothing)
    Dirichlet,  // Prescribed tangential E-field time derivative
    Absorbing,  // Sommerfeld first-order ABC
    PML,        // Perfectly Matched Layer (TODO)
    Port        // Waveguide port condition (TODO)
};

/**
 * @brief Boundary condition manager for Maxwell solver
 * 
 * Manages boundary conditions on different surfaces:
 * - Natural BC: n × (∇×E) = 0 (implicitly satisfied)
 * - Dirichlet BC: ∂E/∂t|_boundary = f(x,t)
 * - Absorbing BC: n × (∇×E) = -√(ε/μ) ∂E/∂t|_boundary
 */
class MaxwellBoundaryConditions
{
public:
    /**
     * @brief Constructor
     * @param abcMarkers Array of boundary attribute markers for ABC
     * @param dbcMarkers Array of boundary attribute markers for Dirichlet BC
     * @param dEdtFunc Time derivative of E-field on Dirichlet boundaries
     * @param etaInvCoef Impedance coefficient √(ε/μ) for ABC
     * @param fespace Finite element space (H(curl))
     */
#ifdef MFEM_USE_MPI
    MaxwellBoundaryConditions(const mfem::Array<int>& abcMarkers,
                              const mfem::Array<int>& dbcMarkers,
                              void (*dEdtFunc)(const mfem::Vector&, double, mfem::Vector&),
                              mfem::Coefficient* etaInvCoef,
                              mfem::ParFiniteElementSpace* fespace);
#else
    MaxwellBoundaryConditions(const mfem::Array<int>& abcMarkers,
                              const mfem::Array<int>& dbcMarkers,
                              void (*dEdtFunc)(const mfem::Vector&, double, mfem::Vector&),
                              mfem::Coefficient* etaInvCoef,
                              mfem::FiniteElementSpace* fespace);
#endif

    /**
     * @brief Destructor
     */
    ~MaxwellBoundaryConditions();

    /**
     * @brief Check if absorbing boundaries are present
     */
    bool hasAbsorbingBC() const { return hasABC_; }

    /**
     * @brief Check if Dirichlet boundaries are present
     */
    bool hasDirichletBC() const { return hasDBC_; }

    /**
     * @brief Get Dirichlet DOF list
     */
    const mfem::Array<int>& getDirichletDofs() const { return dbcDofs_; }

    /**
     * @brief Get ABC bilinear form (for adding to system matrix)
     */
#ifdef MFEM_USE_MPI
    mfem::ParBilinearForm* getAbsorbingForm() const { return abcForm_; }
#else
    mfem::BilinearForm* getAbsorbingForm() const { return abcForm_; }
#endif

    /**
     * @brief Update Dirichlet BC for current time
     * @param time Current simulation time
     */
    void updateDirichletBC(double time);

    /**
     * @brief Apply Dirichlet BC to vector
     * @param vec Vector to apply BC to
     */
    void applyDirichletBC(mfem::Vector& vec) const;

private:
    bool hasABC_;
    bool hasDBC_;

    // Boundary attribute markers
    mfem::Array<int> abcMarker_;
    mfem::Array<int> dbcMarker_;

    // Dirichlet DOFs
    mfem::Array<int> dbcDofs_;

    // Dirichlet BC function
    void (*dEdtFunc_)(const mfem::Vector&, double, mfem::Vector&);

    // Current time for BC evaluation
    double currentTime_;

#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* fespace_;
    mfem::ParBilinearForm* abcForm_;
#else
    mfem::FiniteElementSpace* fespace_;
    mfem::BilinearForm* abcForm_;
#endif

    mfem::VectorCoefficient* dEdtCoef_;
};

} // namespace hpcfem

#endif // HPCFEM_MAXWELL_BOUNDARY_HPP
