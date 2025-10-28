/**
 * @file physics_joule_heating.hpp
 * @brief Coupled Joule heating physics implementation
 * 
 * Implements the coupled electro-thermal system:
 * Block (0,0): -∇·(σ(T)∇V) = 0  (electrostatics)
 * Block (1,0): Q = σ|∇V|²        (Joule heating coupling)
 * Block (1,1): -∇·(κ∇T) = Q     (thermal diffusion)
 * 
 * The system is assembled as a 2×2 BlockOperator:
 * [K_e(T)    0  ] [V]   [0]
 * [C(V)    K_t  ] [T] = [0]
 * 
 * Phase: 3, Step: 3.1
 */

#ifndef HPCFEM_PHYSICS_JOULE_HEATING_HPP
#define HPCFEM_PHYSICS_JOULE_HEATING_HPP

#include "mfem.hpp"
#include <memory>

namespace hpcfem
{

/**
 * @class JouleHeatingPhysics
 * @brief Coupled electro-thermal physics for Joule heating
 * 
 * This class assembles a monolithic 2×2 block system for coupled
 * electrostatics and thermal diffusion with Joule heating coupling.
 * 
 * NOTE: This class does NOT implement PhysicsInterface as the coupled
 * system requires BlockOperator/BlockVector which don't fit the
 * simple interface. This is a specialized coupled physics class.
 */
class JouleHeatingPhysics
{
public:
#ifdef MFEM_USE_MPI
    /**
     * @brief Constructor for parallel execution
     * @param pmesh Parallel mesh (ownership not transferred)
     * @param polynomialOrder Polynomial order for H1 finite elements
     * @param electricalConductivity Electrical conductivity σ [S/m]
     * @param thermalConductivity Thermal conductivity κ [W/(m·K)]
     */
    JouleHeatingPhysics(mfem::ParMesh* pmesh,
                        int polynomialOrder,
                        double electricalConductivity,
                        double thermalConductivity);
#else
    /**
     * @brief Constructor for serial execution
     * @param mesh Mesh (ownership not transferred)
     * @param polynomialOrder Polynomial order for H1 finite elements
     * @param electricalConductivity Electrical conductivity σ [S/m]
     * @param thermalConductivity Thermal conductivity κ [W/(m·K)]
     */
    JouleHeatingPhysics(mfem::Mesh* mesh,
                        int polynomialOrder,
                        double electricalConductivity,
                        double thermalConductivity);
#endif

    /**
     * @brief Destructor - cleans up owned resources
     */
    ~JouleHeatingPhysics();

    /**
     * @brief Assemble the coupled 2×2 block system
     * 
     * Assembles the block operator and block RHS:
     * A = [K_e   0 ]    b = [0]
     *     [C    K_t]        [0]
     * 
     * @param blockOperator Output block operator (must be pre-allocated)
     * @param blockRHS Output block right-hand side
     * @param blockSolution Solution vector for boundary conditions
     * @param essTdofElectric Essential DOFs for electric field
     * @param essTdofThermal Essential DOFs for temperature field
     */
#ifdef MFEM_USE_MPI
    void assemble(mfem::BlockOperator& blockOperator,
                 mfem::BlockVector& blockRHS,
                 mfem::BlockVector& blockSolution,
                 mfem::Array<int>& essTdofElectric,
                 mfem::Array<int>& essTdofThermal);
#else
    void assemble(mfem::BlockOperator& blockOperator,
                 mfem::BlockVector& blockRHS,
                 mfem::BlockVector& blockSolution,
                 mfem::Array<int>& essTdofElectric,
                 mfem::Array<int>& essTdofThermal);
#endif

    /**
     * @brief Get the electric finite element space
     * @return Pointer to electric FE space (ownership retained)
     */
#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* getElectricSpace();
#else
    mfem::FiniteElementSpace* getElectricSpace();
#endif

    /**
     * @brief Get the thermal finite element space
     * @return Pointer to thermal FE space (ownership retained)
     */
#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* getThermalSpace();
#else
    mfem::FiniteElementSpace* getThermalSpace();
#endif

    /**
     * @brief Get block offsets for solution vector assembly
     * @return Array of block offsets [0, n_electric, n_electric+n_thermal]
     */
    const mfem::Array<int>& getBlockOffsets() const;

private:
    // Finite element collections (we own these)
    mfem::H1_FECollection* fecElectric_;
    mfem::H1_FECollection* fecThermal_;
    
    // Finite element spaces (we own these)
#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* fespaceElectric_;
    mfem::ParFiniteElementSpace* fespaceThermal_;
    
    // Bilinear forms for each block
    mfem::ParBilinearForm* bilinearElectric_;  // Block (0,0): K_e
    mfem::ParBilinearForm* bilinearThermal_;   // Block (1,1): K_t
    
    // Linear forms
    mfem::ParLinearForm* linearElectric_;
    mfem::ParLinearForm* linearThermal_;
    
    // Assembled matrices (we own these after ParallelAssemble)
    mfem::HypreParMatrix* matrixK_e_;
    mfem::HypreParMatrix* matrixK_t_;
    mfem::HypreParMatrix* matrixC_;  // Coupling matrix (will be computed later)
#else
    mfem::FiniteElementSpace* fespaceElectric_;
    mfem::FiniteElementSpace* fespaceThermal_;
    
    mfem::BilinearForm* bilinearElectric_;
    mfem::BilinearForm* bilinearThermal_;
    
    mfem::LinearForm* linearElectric_;
    mfem::LinearForm* linearThermal_;
    
    mfem::SparseMatrix* matrixK_e_;
    mfem::SparseMatrix* matrixK_t_;
    mfem::SparseMatrix* matrixC_;
#endif

    // Material properties
    double sigma_;  // Electrical conductivity
    double kappa_;  // Thermal conductivity
    
    // Coefficients (must be members to avoid dangling pointers)
    mfem::ConstantCoefficient* sigmaCoeff_;
    mfem::ConstantCoefficient* kappaCoeff_;
    
    // Block structure
    mfem::Array<int> blockOffsets_;
};

} // namespace hpcfem

#endif // HPCFEM_PHYSICS_JOULE_HEATING_HPP
