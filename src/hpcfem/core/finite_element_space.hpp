#ifndef HPCFEM_CORE_FINITE_ELEMENT_SPACE_HPP
#define HPCFEM_CORE_FINITE_ELEMENT_SPACE_HPP

#include "mfem.hpp"

namespace hpcfem
{

class FiniteElementSpace
{
public:
    enum class Type { H1, L2, ND_VECTOR, RT_VECTOR };

#ifdef MFEM_USE_MPI
    FiniteElementSpace(mfem::ParMesh* pmesh, int order, int dim = 1, Type type = Type::H1);
    mfem::ParFiniteElementSpace* getMfemSpace() const;
#else
    FiniteElementSpace(mfem::Mesh* mesh, int order, int dim = 1, Type type = Type::H1);
    mfem::FiniteElementSpace* getMfemSpace() const;
#endif

    ~FiniteElementSpace();

    mfem::FiniteElementCollection* getFeCollection() const;

private:
    mfem::FiniteElementCollection* fec_ = nullptr; // Owned
#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* fes_ = nullptr; // Owned
    mfem::ParMesh* pmesh_ = nullptr; // Not owned
#else
    mfem::FiniteElementSpace* fes_ = nullptr; // Owned
    mfem::Mesh* mesh_ = nullptr; // Not owned
#endif
    bool isOwner_ = true; // Manages fec_ and fes_
};

} // namespace hpcfem

#endif // HPCFEM_CORE_FINITE_ELEMENT_SPACE_HPP
