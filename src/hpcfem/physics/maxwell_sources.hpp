/**
 * @file maxwell_sources.hpp
 * @brief Current density source term management for Maxwell solver
 * 
 * Handles time-dependent volumetric current density sources J(x,t)
 * for the Maxwell equations.
 */

#ifndef HPCFEM_MAXWELL_SOURCES_HPP
#define HPCFEM_MAXWELL_SOURCES_HPP

#include "mfem.hpp"

namespace hpcfem
{

/**
 * @brief Current density source term manager
 * 
 * Manages time-dependent volumetric current density J(x,t) which
 * appears on the RHS of Ampere's law:
 *   ε ∂E/∂t = ∇×(μ⁻¹B) - σE - J
 */
class MaxwellSourceTerms
{
public:
    /**
     * @brief Constructor with current density function
     * @param jFunc Current density function J(x,t,j_vec)
     * @param fespace Finite element space for current density (H(curl))
     * 
     * If jFunc is nullptr, no source term is present.
     */
#ifdef MFEM_USE_MPI
    MaxwellSourceTerms(void (*jFunc)(const mfem::Vector&, double, mfem::Vector&),
                       mfem::ParFiniteElementSpace* fespace);
#else
    MaxwellSourceTerms(void (*jFunc)(const mfem::Vector&, double, mfem::Vector&),
                       mfem::FiniteElementSpace* fespace);
#endif

    /**
     * @brief Destructor
     */
    ~MaxwellSourceTerms();

    /**
     * @brief Check if source term exists
     * @return True if J(x,t) is present
     */
    bool hasSource() const { return jFunc_ != nullptr; }

    /**
     * @brief Update source term for current time
     * @param time Current simulation time
     */
    void updateTime(double time);

    /**
     * @brief Get the dual of current density (for RHS assembly)
     * @return Pointer to linear form containing -(J, F) where F are test functions
     */
#ifdef MFEM_USE_MPI
    mfem::ParLinearForm* getDualForm() const { return jDual_; }
#else
    mfem::LinearForm* getDualForm() const { return jDual_; }
#endif

    /**
     * @brief Get the current density as a grid function
     * @return Pointer to grid function representing J(x,t)
     */
#ifdef MFEM_USE_MPI
    mfem::ParGridFunction* getCurrentDensity() const { return jField_; }
#else
    mfem::GridFunction* getCurrentDensity() const { return jField_; }
#endif

private:
    // Current density function
    void (*jFunc_)(const mfem::Vector&, double, mfem::Vector&);

    // Current time
    double currentTime_;

    // Finite element space
#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* fespace_;
    mfem::ParLinearForm* jDual_;
    mfem::ParGridFunction* jField_;
#else
    mfem::FiniteElementSpace* fespace_;
    mfem::LinearForm* jDual_;
    mfem::GridFunction* jField_;
#endif

    // Vector coefficient for J(x,t)
    mfem::VectorCoefficient* jCoef_;
};

} // namespace hpcfem

#endif // HPCFEM_MAXWELL_SOURCES_HPP
