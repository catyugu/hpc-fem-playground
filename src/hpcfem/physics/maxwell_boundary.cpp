/**
 * @file maxwell_boundary.cpp
 * @brief Implementation of MaxwellBoundaryConditions class
 */

#include "maxwell_boundary.hpp"

namespace hpcfem
{

// Helper class for time-dependent boundary condition
class TimeDependentBCCoefficient : public mfem::VectorCoefficient
{
public:
    TimeDependentBCCoefficient(
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

MaxwellBoundaryConditions::MaxwellBoundaryConditions(
    const mfem::Array<int>& abcMarkers,
    const mfem::Array<int>& dbcMarkers,
    void (*dEdtFunc)(const mfem::Vector&, double, mfem::Vector&),
    mfem::Coefficient* etaInvCoef,
    mfem::ParFiniteElementSpace* fespace)
    : hasABC_(abcMarkers.Size() > 0),
      hasDBC_(dbcMarkers.Size() > 0),
      dEdtFunc_(dEdtFunc),
      currentTime_(0.0),
      fespace_(fespace),
      abcForm_(nullptr),
      dEdtCoef_(nullptr)
{
    // Setup ABC marker array
    if (hasABC_)
    {
        int numBdrAttr = fespace_->GetMesh()->bdr_attributes.Max();
        abcMarker_.SetSize(numBdrAttr);
        abcMarker_ = 0;
        for (int i = 0; i < abcMarkers.Size(); i++)
        {
            abcMarker_[abcMarkers[i] - 1] = 1;
        }

        // Create ABC bilinear form
        abcForm_ = new mfem::ParBilinearForm(fespace_);
        abcForm_->AddBoundaryIntegrator(
            new mfem::VectorFEMassIntegrator(*etaInvCoef), abcMarker_);
        abcForm_->Assemble();
        abcForm_->Finalize();
    }

    // Setup Dirichlet BC marker and DOFs
    if (hasDBC_)
    {
        int numBdrAttr = fespace_->GetMesh()->bdr_attributes.Max();
        dbcMarker_.SetSize(numBdrAttr);
        dbcMarker_ = 0;
        for (int i = 0; i < dbcMarkers.Size(); i++)
        {
            dbcMarker_[dbcMarkers[i] - 1] = 1;
        }

        // Get essential DOFs
        fespace_->GetEssentialTrueDofs(dbcMarker_, dbcDofs_);

        // Create BC coefficient if function provided
        if (dEdtFunc_ != nullptr)
        {
            int dim = fespace_->GetMesh()->SpaceDimension();
            dEdtCoef_ = new TimeDependentBCCoefficient(dim, dEdtFunc_);
        }
    }
}

#else

MaxwellBoundaryConditions::MaxwellBoundaryConditions(
    const mfem::Array<int>& abcMarkers,
    const mfem::Array<int>& dbcMarkers,
    void (*dEdtFunc)(const mfem::Vector&, double, mfem::Vector&),
    mfem::Coefficient* etaInvCoef,
    mfem::FiniteElementSpace* fespace)
    : hasABC_(abcMarkers.Size() > 0),
      hasDBC_(dbcMarkers.Size() > 0),
      dEdtFunc_(dEdtFunc),
      currentTime_(0.0),
      fespace_(fespace),
      abcForm_(nullptr),
      dEdtCoef_(nullptr)
{
    // Setup ABC marker array
    if (hasABC_)
    {
        int numBdrAttr = fespace_->GetMesh()->bdr_attributes.Max();
        abcMarker_.SetSize(numBdrAttr);
        abcMarker_ = 0;
        for (int i = 0; i < abcMarkers.Size(); i++)
        {
            abcMarker_[abcMarkers[i] - 1] = 1;
        }

        // Create ABC bilinear form
        abcForm_ = new mfem::BilinearForm(fespace_);
        abcForm_->AddBoundaryIntegrator(
            new mfem::VectorFEMassIntegrator(*etaInvCoef), abcMarker_);
        abcForm_->Assemble();
        abcForm_->Finalize();
    }

    // Setup Dirichlet BC marker and DOFs
    if (hasDBC_)
    {
        int numBdrAttr = fespace_->GetMesh()->bdr_attributes.Max();
        dbcMarker_.SetSize(numBdrAttr);
        dbcMarker_ = 0;
        for (int i = 0; i < dbcMarkers.Size(); i++)
        {
            dbcMarker_[dbcMarkers[i] - 1] = 1;
        }

        // Get essential DOFs
        fespace_->GetEssentialTrueDofs(dbcMarker_, dbcDofs_);

        // Create BC coefficient if function provided
        if (dEdtFunc_ != nullptr)
        {
            int dim = fespace_->GetMesh()->SpaceDimension();
            dEdtCoef_ = new TimeDependentBCCoefficient(dim, dEdtFunc_);
        }
    }
}

#endif

MaxwellBoundaryConditions::~MaxwellBoundaryConditions()
{
    delete abcForm_;
    delete dEdtCoef_;
}

void MaxwellBoundaryConditions::updateDirichletBC(double time)
{
    if (!hasDBC_ || dEdtFunc_ == nullptr) return;

    currentTime_ = time;

    TimeDependentBCCoefficient* tdCoef = 
        dynamic_cast<TimeDependentBCCoefficient*>(dEdtCoef_);
    if (tdCoef != nullptr)
    {
        tdCoef->setTime(time);
    }
}

void MaxwellBoundaryConditions::applyDirichletBC(mfem::Vector& vec) const
{
    if (!hasDBC_) return;

    // Set Dirichlet DOFs to zero (homogeneous BC)
    for (int i = 0; i < dbcDofs_.Size(); i++)
    {
        vec(dbcDofs_[i]) = 0.0;
    }
}

} // namespace hpcfem
