/**
 * @file physics_electrostatics.hpp
 * @brief Electrostatics physics implementation
 * 
 * Implements the electrostatics equation: ∇·(σ∇φ) = f
 * where σ is the electrical conductivity and f is the source term.
 */

#ifndef HPCFEM_PHYSICS_ELECTROSTATICS_HPP
#define HPCFEM_PHYSICS_ELECTROSTATICS_HPP

#include "physics_interface.hpp"
#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class ElectrostaticsPhysics
 * @brief Concrete implementation for electrostatics problems
 * 
 * This class handles the assembly of the electrostatics problem:
 * - Bilinear form: a(u,v) = ∫ σ∇u·∇v dx
 * - Linear form: l(v) = ∫ f·v dx
 */
class ElectrostaticsPhysics : public PhysicsInterface
{
public:
#ifdef MFEM_USE_MPI
    /**
     * @brief Constructor for parallel execution
     * @param pmesh Parallel mesh (ownership not transferred)
     * @param polynomialOrder Polynomial order for H1 finite elements
     * @param conductivity Electrical conductivity coefficient
     * @param sourceCoeff Source term coefficient
     * @param boundaryCoeff Boundary condition coefficient
     */
    ElectrostaticsPhysics(mfem::ParMesh* pmesh,
                         int polynomialOrder,
                         mfem::Coefficient* conductivity,
                         mfem::Coefficient* sourceCoeff,
                         mfem::Coefficient* boundaryCoeff);
#else
    /**
     * @brief Constructor for serial execution
     * @param mesh Mesh (ownership not transferred)
     * @param polynomialOrder Polynomial order for H1 finite elements
     * @param conductivity Electrical conductivity coefficient
     * @param sourceCoeff Source term coefficient
     * @param boundaryCoeff Boundary condition coefficient
     */
    ElectrostaticsPhysics(mfem::Mesh* mesh,
                         int polynomialOrder,
                         mfem::Coefficient* conductivity,
                         mfem::Coefficient* sourceCoeff,
                         mfem::Coefficient* boundaryCoeff);
#endif

    /**
     * @brief Destructor - cleans up owned resources
     */
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
