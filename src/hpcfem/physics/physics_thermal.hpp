/**
 * @file physics/physics_thermal.hpp
 * @brief Thermal diffusion physics implementation
 */

#ifndef HPCFEM_PHYSICS_THERMAL_HPP
#define HPCFEM_PHYSICS_THERMAL_HPP

#include "hpcfem/core/physics_interface.hpp"
#include "mfem.hpp"

namespace hpcfem
{

class ThermalPhysics : public PhysicsInterface
{
public:
#ifdef MFEM_USE_MPI
    ThermalPhysics(mfem::ParMesh* pmesh,
                   int polynomialOrder,
                   mfem::Coefficient* thermalConductivity,
                   mfem::Coefficient* heatSource,
                   mfem::Coefficient* boundaryCoeff);
#else
    ThermalPhysics(mfem::Mesh* mesh,
                   int polynomialOrder,
                   mfem::Coefficient* thermalConductivity,
                   mfem::Coefficient* heatSource,
                   mfem::Coefficient* boundaryCoeff);
#endif

    ~ThermalPhysics() override;

#ifdef MFEM_USE_MPI
    void assemble(mfem::HypreParMatrix& A,
                 mfem::Vector& b,
                 mfem::Vector& x,
                 mfem::Array<int>& essTdofList) override;

    mfem::ParFiniteElementSpace* getFiniteElementSpace() override;
#else
    void assemble(mfem::SparseMatrix& A,
                 mfem::Vector& b,
                 mfem::Vector& x,
                 mfem::Array<int>& essTdofList) override;

    mfem::FiniteElementSpace* getFiniteElementSpace() override;
#endif

private:
    mfem::H1_FECollection* fec_;
#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* fespace_;
    mfem::ParBilinearForm* bilinearForm_;
    mfem::ParLinearForm* linearForm_;
#else
    mfem::FiniteElementSpace* fespace_;
    mfem::BilinearForm* bilinearForm_;
    mfem::LinearForm* linearForm_;
#endif

    mfem::Coefficient* thermalConductivity_;
    mfem::Coefficient* heatSource_;
    mfem::Coefficient* boundaryCoeff_;
};

} // namespace hpcfem

#endif // HPCFEM_PHYSICS_THERMAL_HPP
