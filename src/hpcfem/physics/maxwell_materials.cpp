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
    if (epsilonFunc_ != nullptr || muInvFunc_ != nullptr)
    {
        // For spatially varying materials, use product and transform
        // Create product coefficient: ε·μ⁻¹
        mfem::ProductCoefficient* epsTimesMusInv = 
            new mfem::ProductCoefficient(*epsilonCoef_, *muInvCoef_);
        
        // Apply sqrt transformation: √(ε·μ⁻¹) = η⁻¹
        etaInvCoef_ = new mfem::TransformedCoefficient(
            epsTimesMusInv, [](double x) { return sqrt(x); });
        
        // Note: TransformedCoefficient takes ownership of the coefficient
    }
    else
    {
        // For uniform materials, can compute directly
        double eps = EPSILON_0;
        double muInv = 1.0 / MU_0;
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
