/**
 * @file maxwell_materials.cpp
 * @brief Implementation of MaxwellMaterialProperties class
 */

#include "maxwell_materials.hpp"
#include <cmath>

namespace hpcfem
{

// Helper function to compute sqrt(epsilon/mu) for impedance
double computeImpedanceInv(double (*epsilonFunc)(const mfem::Vector&),
                           double (*muInvFunc)(const mfem::Vector&),
                           const mfem::Vector& x)
{
    double eps = epsilonFunc(x);
    double muInv = muInvFunc(x);
    return sqrt(eps * muInv);
}

MaxwellMaterialProperties::MaxwellMaterialProperties(
    double (*epsilonFunc)(const mfem::Vector&),
    double (*muInvFunc)(const mfem::Vector&),
    double (*sigmaFunc)(const mfem::Vector&))
    : epsilonFunc_(epsilonFunc),
      muInvFunc_(muInvFunc),
      sigmaFunc_(sigmaFunc),
      epsilonCoef_(nullptr),
      muInvCoef_(nullptr),
      sigmaCoef_(nullptr),
      etaInvCoef_(nullptr)
{
    // Create permittivity coefficient
    if (epsilonFunc_ == nullptr)
    {
        // Default to free space
        epsilonCoef_ = new mfem::ConstantCoefficient(EPSILON_0);
    }
    else
    {
        epsilonCoef_ = new mfem::FunctionCoefficient(epsilonFunc_);
    }

    // Create inverse permeability coefficient
    if (muInvFunc_ == nullptr)
    {
        // Default to free space
        muInvCoef_ = new mfem::ConstantCoefficient(1.0 / MU_0);
    }
    else
    {
        muInvCoef_ = new mfem::FunctionCoefficient(muInvFunc_);
    }

    // Create conductivity coefficient (only if lossy)
    if (sigmaFunc_ != nullptr)
    {
        sigmaCoef_ = new mfem::FunctionCoefficient(sigmaFunc_);
    }

    // Create impedance coefficient for ABC
    // η⁻¹ = √(ε/μ) = √(ε·μ⁻¹)
    if (epsilonFunc_ != nullptr && muInvFunc_ != nullptr)
    {
        // For spatially varying materials, need to compute impedance
        // This is a simplified implementation - proper implementation
        // would use TransformedCoefficient
        // TODO: Implement proper spatially-varying impedance
        etaInvCoef_ = nullptr;
    }
    else
    {
        // For uniform materials, can compute directly
        double eps = (epsilonFunc_ == nullptr) ? EPSILON_0 : 
                     epsilonFunc_(mfem::Vector(3));
        double muInv = (muInvFunc_ == nullptr) ? (1.0 / MU_0) : 
                       muInvFunc_(mfem::Vector(3));
        double etaInv = sqrt(eps * muInv);
        etaInvCoef_ = new mfem::ConstantCoefficient(etaInv);
    }
}

MaxwellMaterialProperties::~MaxwellMaterialProperties()
{
    delete epsilonCoef_;
    delete muInvCoef_;
    delete sigmaCoef_;
    delete etaInvCoef_;
}

} // namespace hpcfem
