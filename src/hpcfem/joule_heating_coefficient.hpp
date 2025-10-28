/**
 * @file joule_heating_coefficient.hpp
 * @brief Coefficient for Joule heating source term Q = σ|∇V|²
 * 
 * Implements a grid function coefficient that computes the Joule heating
 * power density from the electric field gradient.
 * 
 * Phase: 3, Step: 3.2 (nonlinear coupling)
 */

#ifndef HPCFEM_JOULE_HEATING_COEFFICIENT_HPP
#define HPCFEM_JOULE_HEATING_COEFFICIENT_HPP

#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class JouleHeatingCoefficient
 * @brief Coefficient that evaluates Q = σ|∇V|²
 * 
 * This coefficient wraps a GridFunction for voltage V and computes
 * the Joule heating power density at each quadrature point.
 * 
 * The formula is:
 *   Q(x) = σ * |∇V(x)|²
 * 
 * where:
 *   - σ is the electrical conductivity [S/m]
 *   - ∇V is the electric field gradient [V/m]
 *   - Q is the volumetric heating rate [W/m³]
 */
#ifdef MFEM_USE_MPI
class JouleHeatingCoefficient : public mfem::Coefficient
{
public:
    /**
     * @brief Constructor
     * @param voltage Grid function for voltage V (not owned)
     * @param sigma Electrical conductivity [S/m]
     */
    JouleHeatingCoefficient(mfem::ParGridFunction* voltage, double sigma)
        : voltage_(voltage), sigma_(sigma) {}

    /**
     * @brief Evaluate Q = σ|∇V|² at a quadrature point
     * @param T Element transformation
     * @param ip Integration point
     * @return Joule heating power density [W/m³]
     */
    virtual double Eval(mfem::ElementTransformation& T,
                       const mfem::IntegrationPoint& ip) override
    {
        // Get gradient of voltage at this point
        mfem::Vector grad_v;
        voltage_->GetGradient(T, grad_v);
        
        // Compute |∇V|²
        double grad_v_squared = grad_v * grad_v;
        
        // Return Q = σ|∇V|²
        return sigma_ * grad_v_squared;
    }

private:
    mfem::ParGridFunction* voltage_;  ///< Voltage field V (not owned)
    double sigma_;                    ///< Electrical conductivity σ
};
#else
class JouleHeatingCoefficient : public mfem::Coefficient
{
public:
    JouleHeatingCoefficient(mfem::GridFunction* voltage, double sigma)
        : voltage_(voltage), sigma_(sigma) {}

    virtual double Eval(mfem::ElementTransformation& T,
                       const mfem::IntegrationPoint& ip) override
    {
        mfem::Vector grad_v;
        voltage_->GetGradient(T, grad_v);
        double grad_v_squared = grad_v * grad_v;
        return sigma_ * grad_v_squared;
    }

private:
    mfem::GridFunction* voltage_;
    double sigma_;
};
#endif

} // namespace hpcfem

#endif // HPCFEM_JOULE_HEATING_COEFFICIENT_HPP
