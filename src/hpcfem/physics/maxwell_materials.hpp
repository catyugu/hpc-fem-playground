/**
 * @file maxwell_materials.hpp
 * @brief Material property management for Maxwell time-domain solver
 * 
 * Handles electric permittivity, magnetic permeability, and conductivity
 * coefficients for the Maxwell equations. Supports both scalar and 
 * spatially-varying material properties.
 */

#ifndef HPCFEM_MAXWELL_MATERIALS_HPP
#define HPCFEM_MAXWELL_MATERIALS_HPP

#include "mfem.hpp"

namespace hpcfem
{

/**
 * @brief Material property manager for Maxwell solver
 * 
 * Manages three key material properties:
 * - Electric permittivity ε (F/m)
 * - Magnetic permeability μ (H/m), stored as inverse μ⁻¹
 * - Electric conductivity σ (S/m)
 * 
 * All properties can be spatially varying through function pointers.
 */
class MaxwellMaterialProperties
{
public:
    /**
     * @brief Constructor with material property functions
     * @param epsilonFunc Electric permittivity function ε(x)
     * @param muInvFunc Inverse magnetic permeability function μ⁻¹(x)
     * @param sigmaFunc Electric conductivity function σ(x), can be nullptr
     * 
     * If functions are nullptr, defaults to free space values:
     * - ε₀ = 8.8541878176e-12 F/m
     * - μ₀⁻¹ = 1/(4π×10⁻⁷) H⁻¹m⁻¹
     * - σ = 0 S/m (lossless)
     */
    MaxwellMaterialProperties(double (*epsilonFunc)(const mfem::Vector&),
                             double (*muInvFunc)(const mfem::Vector&),
                             double (*sigmaFunc)(const mfem::Vector&));

    /**
     * @brief Destructor
     */
    ~MaxwellMaterialProperties();

    /**
     * @brief Get the permittivity coefficient
     * @return Pointer to MFEM Coefficient for ε(x)
     */
    mfem::Coefficient* getPermittivityCoefficient() const { return epsilonCoef_; }

    /**
     * @brief Get the inverse permeability coefficient
     * @return Pointer to MFEM Coefficient for μ⁻¹(x)
     */
    mfem::Coefficient* getPermeabilityInvCoefficient() const { return muInvCoef_; }

    /**
     * @brief Get the conductivity coefficient
     * @return Pointer to MFEM Coefficient for σ(x), or nullptr if lossless
     */
    mfem::Coefficient* getConductivityCoefficient() const { return sigmaCoef_; }

    /**
     * @brief Check if material is lossy (has conductivity)
     * @return True if σ > 0 anywhere
     */
    bool isLossy() const { return sigmaCoef_ != nullptr; }

    /**
     * @brief Get the impedance coefficient for ABC
     * @return Pointer to MFEM Coefficient for η⁻¹ = √(ε/μ), or nullptr
     */
    mfem::Coefficient* getImpedanceInvCoefficient() const { return etaInvCoef_; }

private:
    // Material property functions
    double (*epsilonFunc_)(const mfem::Vector&);
    double (*muInvFunc_)(const mfem::Vector&);
    double (*sigmaFunc_)(const mfem::Vector&);

    // MFEM coefficients
    mfem::Coefficient* epsilonCoef_;
    mfem::Coefficient* muInvCoef_;
    mfem::Coefficient* sigmaCoef_;
    mfem::Coefficient* etaInvCoef_;

    // Free space constants
    static constexpr double EPSILON_0 = 8.8541878176e-12;  // F/m
    static constexpr double MU_0 = 4.0e-7 * M_PI;           // H/m
};

} // namespace hpcfem

#endif // HPCFEM_MAXWELL_MATERIALS_HPP
