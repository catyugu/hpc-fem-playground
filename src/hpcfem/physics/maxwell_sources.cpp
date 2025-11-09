/**
 * @file maxwell_sources.cpp
 * @brief Implementation of MaxwellSourceTerms class
 */

#include "maxwell_sources.hpp"

namespace hpcfem
{

// Helper class for time-dependent vector coefficient
class TimeDependentVectorCoefficient : public mfem::VectorCoefficient
{
public:
    TimeDependentVectorCoefficient(
        int dim,
        void (*func)(const mfem::Vector&, double, mfem::Vector&))
        : mfem::VectorCoefficient(dim), func_(func), time_(0.0) {}

    void setTime(double t) { time_ = t; }

    virtual void Eval(mfem::Vector& V, mfem::ElementTransformation& T,
                     const mfem::IntegrationPoint& ip)
    {
        double x[3];
        mfem::Vector transip(x, 3);
        T.Transform(ip, transip);
        func_(transip, time_, V);
    }

private:
    void (*func_)(const mfem::Vector&, double, mfem::Vector&);
    double time_;
};

#ifdef MFEM_USE_MPI

MaxwellSourceTerms::MaxwellSourceTerms(
    void (*jFunc)(const mfem::Vector&, double, mfem::Vector&),
    mfem::ParFiniteElementSpace* fespace)
    : jFunc_(jFunc),
      currentTime_(0.0),
      fespace_(fespace),
      jDual_(nullptr),
      jField_(nullptr),
      jCoef_(nullptr)
{
    if (jFunc_ != nullptr)
    {
        // Create vector coefficient
        int dim = fespace_->GetMesh()->SpaceDimension();
        jCoef_ = new TimeDependentVectorCoefficient(dim, jFunc_);

        // Create grid function for visualization
        jField_ = new mfem::ParGridFunction(fespace_);

        // Create linear form for RHS assembly
        jDual_ = new mfem::ParLinearForm(fespace_);
        jDual_->AddDomainIntegrator(
            new mfem::VectorFEDomainLFIntegrator(*jCoef_));
    }
}

#else

MaxwellSourceTerms::MaxwellSourceTerms(
    void (*jFunc)(const mfem::Vector&, double, mfem::Vector&),
    mfem::FiniteElementSpace* fespace)
    : jFunc_(jFunc),
      currentTime_(0.0),
      fespace_(fespace),
      jDual_(nullptr),
      jField_(nullptr),
      jCoef_(nullptr)
{
    if (jFunc_ != nullptr)
    {
        // Create vector coefficient
        int dim = fespace_->GetMesh()->SpaceDimension();
        jCoef_ = new TimeDependentVectorCoefficient(dim, jFunc_);

        // Create grid function for visualization
        jField_ = new mfem::GridFunction(fespace_);

        // Create linear form for RHS assembly
        jDual_ = new mfem::LinearForm(fespace_);
        jDual_->AddDomainIntegrator(
            new mfem::VectorFEDomainLFIntegrator(*jCoef_));
    }
}

#endif

MaxwellSourceTerms::~MaxwellSourceTerms()
{
    delete jCoef_;
    delete jField_;
    delete jDual_;
}

void MaxwellSourceTerms::updateTime(double time)
{
    if (jFunc_ == nullptr) return;

    currentTime_ = time;

    // Update coefficient time
    TimeDependentVectorCoefficient* tdCoef = 
        dynamic_cast<TimeDependentVectorCoefficient*>(jCoef_);
    if (tdCoef != nullptr)
    {
        tdCoef->setTime(time);
    }

    // Project onto grid function for visualization
    jField_->ProjectCoefficient(*jCoef_);

    // Reassemble linear form
    jDual_->Assemble();
}

} // namespace hpcfem
