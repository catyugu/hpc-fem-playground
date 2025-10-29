#include "hpcfem/core/finite_element_space.hpp"
#include <stdexcept>

namespace hpcfem
{

#ifdef MFEM_USE_MPI
FiniteElementSpace::FiniteElementSpace(mfem::ParMesh* pmesh, int order, int dim, Type type)
    : pmesh_(pmesh)
{
    if (!pmesh_) { throw std::runtime_error("ParMesh pointer is null"); }

    switch (type)
    {
    case Type::H1:
        fec_ = new mfem::H1_FECollection(order, dim);
        break;
    case Type::L2:
        fec_ = new mfem::L2_FECollection(order, dim);
        break;
    case Type::ND_VECTOR:
        fec_ = new mfem::ND_FECollection(order, dim);
        break;
    case Type::RT_VECTOR:
        fec_ = new mfem::RT_FECollection(order, dim);
        break;
    default:
        fec_ = new mfem::H1_FECollection(order, dim);
    }

    fes_ = new mfem::ParFiniteElementSpace(pmesh_, fec_, dim);
}

mfem::ParFiniteElementSpace* FiniteElementSpace::getMfemSpace() const { return fes_; }

#else
FiniteElementSpace::FiniteElementSpace(mfem::Mesh* mesh, int order, int dim, Type type)
    : mesh_(mesh)
{
    if (!mesh_) { throw std::runtime_error("Mesh pointer is null"); }

    switch (type)
    {
    case Type::H1:
        fec_ = new mfem::H1_FECollection(order, dim);
        break;
    case Type::L2:
        fec_ = new mfem::L2_FECollection(order, dim);
        break;
    case Type::ND_VECTOR:
        fec_ = new mfem::ND_FECollection(order, dim);
        break;
    case Type::RT_VECTOR:
        fec_ = new mfem::RT_FECollection(order, dim);
        break;
    default:
        fec_ = new mfem::H1_FECollection(order, dim);
    }

    fes_ = new mfem::FiniteElementSpace(mesh_, fec_, dim);
}

mfem::FiniteElementSpace* FiniteElementSpace::getMfemSpace() const { return fes_; }

#endif

FiniteElementSpace::~FiniteElementSpace()
{
#ifdef MFEM_USE_MPI
    delete fes_;
    delete fec_;
    fes_ = nullptr;
    fec_ = nullptr;
#else
    delete fes_;
    delete fec_;
    fes_ = nullptr;
    fec_ = nullptr;
#endif
}

mfem::FiniteElementCollection* FiniteElementSpace::getFeCollection() const { return fec_; }

} // namespace hpcfem
