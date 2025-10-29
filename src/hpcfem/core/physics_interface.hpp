/**
 * @file core/physics_interface.hpp
 * @brief Abstract interface for physics modules in the FEM framework
 */

#ifndef HPCFEM_CORE_PHYSICS_INTERFACE_HPP
#define HPCFEM_CORE_PHYSICS_INTERFACE_HPP

#include "mfem.hpp"

namespace hpcfem
{

class PhysicsInterface
{
public:
    virtual ~PhysicsInterface() = default;

#ifdef MFEM_USE_MPI
    virtual void assemble(mfem::HypreParMatrix& A,
                         mfem::Vector& b,
                         mfem::Vector& x,
                         mfem::Array<int>& essTdofList) = 0;

    virtual mfem::ParFiniteElementSpace* getFiniteElementSpace() = 0;
#else
    virtual void assemble(mfem::SparseMatrix& A,
                         mfem::Vector& b,
                         mfem::Vector& x,
                         mfem::Array<int>& essTdofList) = 0;

    virtual mfem::FiniteElementSpace* getFiniteElementSpace() = 0;
#endif
};

} // namespace hpcfem

#endif // HPCFEM_CORE_PHYSICS_INTERFACE_HPP
