/**
 * @file physics_interface.hpp
 * @brief Abstract interface for physics modules in the FEM framework
 * 
 * This interface defines the contract for all physics implementations.
 * Each physics module is responsible for assembling the bilinear and linear
 * forms for its governing equations.
 */

#ifndef HPCFEM_PHYSICS_INTERFACE_HPP
#define HPCFEM_PHYSICS_INTERFACE_HPP

#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class PhysicsInterface
 * @brief Abstract base class for all physics modules
 * 
 * This class defines the interface that all physics implementations must follow.
 * Physics modules are responsible for:
 * - Setting up finite element spaces
 * - Assembling bilinear forms (system matrix)
 * - Assembling linear forms (right-hand side)
 * - Applying boundary conditions
 */
class PhysicsInterface
{
public:
    /**
     * @brief Virtual destructor
     */
    virtual ~PhysicsInterface() = default;

#ifdef MFEM_USE_MPI
    /**
     * @brief Assemble the linear system for parallel execution
     * @param A Output system matrix (parallel)
     * @param b Output right-hand side vector
     * @param x Solution vector (used for boundary conditions)
     * @param essTdofList List of essential (Dirichlet) true DOFs
     */
    virtual void assemble(mfem::HypreParMatrix& A,
                         mfem::Vector& b,
                         mfem::Vector& x,
                         mfem::Array<int>& essTdofList) = 0;
#else
    /**
     * @brief Assemble the linear system for serial execution
     * @param A Output system matrix (serial)
     * @param b Output right-hand side vector
     * @param x Solution vector (used for boundary conditions)
     * @param essTdofList List of essential (Dirichlet) true DOFs
     */
    virtual void assemble(mfem::SparseMatrix& A,
                         mfem::Vector& b,
                         mfem::Vector& x,
                         mfem::Array<int>& essTdofList) = 0;
#endif

    /**
     * @brief Get the finite element space for this physics
     * @return Pointer to the finite element space (ownership retained)
     */
#ifdef MFEM_USE_MPI
    virtual mfem::ParFiniteElementSpace* getFiniteElementSpace() = 0;
#else
    virtual mfem::FiniteElementSpace* getFiniteElementSpace() = 0;
#endif
};

} // namespace hpcfem

#endif // HPCFEM_PHYSICS_INTERFACE_HPP
