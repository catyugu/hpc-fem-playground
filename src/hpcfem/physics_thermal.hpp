/**
 * @file physics_thermal.hpp
 * @brief Thermal diffusion physics implementation
 * 
 * Implements the steady-state heat equation: ∇·(κ∇T) = Q
 * where κ is the thermal conductivity and Q is the heat source term.
 * 
 */

#ifndef HPCFEM_PHYSICS_THERMAL_HPP
#define HPCFEM_PHYSICS_THERMAL_HPP

#include "physics_interface.hpp"
#include "mfem.hpp"

namespace hpcfem
{

/**
 * @class ThermalPhysics
 * @brief Concrete implementation for thermal diffusion problems
 * 
 * This class handles the assembly of the steady-state thermal problem:
 * - Bilinear form: a(T,v) = ∫ κ∇T·∇v dx
 * - Linear form: l(v) = ∫ Q·v dx
 */
class ThermalPhysics : public PhysicsInterface
{
public:
#ifdef MFEM_USE_MPI
    /**
     * @brief Constructor for parallel execution
     * @param pmesh Parallel mesh (ownership not transferred)
     * @param polynomialOrder Polynomial order for H1 finite elements
     * @param thermalConductivity Thermal conductivity coefficient κ
     * @param heatSource Heat source coefficient Q
     * @param boundaryCoeff Boundary temperature coefficient
     */
    ThermalPhysics(mfem::ParMesh* pmesh,
                   int polynomialOrder,
                   mfem::Coefficient* thermalConductivity,
                   mfem::Coefficient* heatSource,
                   mfem::Coefficient* boundaryCoeff);
#else
    /**
     * @brief Constructor for serial execution
     * @param mesh Mesh (ownership not transferred)
     * @param polynomialOrder Polynomial order for H1 finite elements
     * @param thermalConductivity Thermal conductivity coefficient κ
     * @param heatSource Heat source coefficient Q
     * @param boundaryCoeff Boundary temperature coefficient
     */
    ThermalPhysics(mfem::Mesh* mesh,
                   int polynomialOrder,
                   mfem::Coefficient* thermalConductivity,
                   mfem::Coefficient* heatSource,
                   mfem::Coefficient* boundaryCoeff);
#endif

    /**
     * @brief Destructor - cleans up owned resources
     */
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
