/**
 * @file physics/physics_electrostatics.hpp
 * @brief Electrostatics physics implementation
 */

#ifndef HPCFEM_PHYSICS_ELECTROSTATICS_HPP
#define HPCFEM_PHYSICS_ELECTROSTATICS_HPP

#include "hpcfem/core/physics_interface.hpp"
#include "mfem.hpp"

namespace hpcfem
{

class ElectrostaticsPhysics : public PhysicsInterface
{
public:
#ifdef MFEM_USE_MPI
    ElectrostaticsPhysics(mfem::ParMesh* pmesh,
                         int polynomialOrder,
                         mfem::Coefficient* conductivity,
                         mfem::Coefficient* sourceCoeff,
                         mfem::Coefficient* boundaryCoeff);
#else
    ElectrostaticsPhysics(mfem::Mesh* mesh,
                         int polynomialOrder,
                         mfem::Coefficient* conductivity,
                         mfem::Coefficient* sourceCoeff,
                         mfem::Coefficient* boundaryCoeff);
#endif

    ~ElectrostaticsPhysics() override;

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

    mfem::Coefficient* conductivity_;
    mfem::Coefficient* sourceCoeff_;
    mfem::Coefficient* boundaryCoeff_;
};

} // namespace hpcfem

#endif // HPCFEM_PHYSICS_ELECTROSTATICS_HPP
