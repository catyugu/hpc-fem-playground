/**
 * @file physics/physics_cavity_eigen.hpp
 * @brief Solves the 3D Maxwell eigenvalue problem for cavity modes
 * 
 * This class implements the 3D electromagnetic eigenvalue solver using
 * Nedelec (edge) elements. It solves: ∇×∇×E = λE with PEC boundary conditions.
 * This is the proper formulation for cavity analysis and mode calculation.
 * 
 * Based on MFEM example ex13p.cpp
 */

#ifndef HPCFEM_PHYSICS_WAVEGUIDE_EIGEN_HPP
#define HPCFEM_PHYSICS_WAVEGUIDE_EIGEN_HPP

#include "mfem.hpp"
#include <vector>

namespace hpcfem
{

/**
 * @brief Solves the 3D Maxwell eigenvalue problem for electromagnetic modes
 * @details This class solves ∇×∇×E = λE using Nedelec finite elements
 * in 3D cavities. The eigenvalues λ represent (k_c)² where k_c is the
 * cutoff wavenumber. This is used for mode analysis and
 * S-parameter port definitions.
 */
class PhysicsCavityEigen
{
public:
#ifdef MFEM_USE_MPI
    /**
     * @brief Constructor for parallel execution
     * @param pmesh Parallel mesh representing the 3D cavity or waveguide section
     * @param order Polynomial order for Nedelec finite element space
     */
    PhysicsCavityEigen(mfem::ParMesh* pmesh, int order);
#else
    /**
     * @brief Constructor for serial execution
     * @param mesh Mesh representing the 3D cavity or waveguide section
     * @param order Polynomial order for Nedelec finite element space
     */
    PhysicsCavityEigen(mfem::Mesh* mesh, int order);
#endif

    /**
     * @brief Destructor
     */
    ~PhysicsCavityEigen();

    /**
     * @brief Solves the generalized eigenvalue problem ∇×∇×E = λ·M·E
     * @param numModes The number of eigenvalues to compute
     * @return A vector of eigenvalues (λ = k_c²)
     * @details Uses HypreAME or HypreLOBPCG with HypreAMS preconditioner
     */
    std::vector<double> solveEigenvalues(int numModes);

    /**
     * @brief Get the computed eigenvectors (mode shapes)
     * @return Vector of eigenvectors representing electric field modes
     */
    std::vector<mfem::Vector> getEigenvectors() const;

    /**
     * @brief Get a specific mode shape as a grid function
     * @param modeIndex Index of the mode (0 = lowest eigenvalue)
     * @return Pointer to grid function representing the E-field mode
     * @details Caller is responsible for deleting the returned pointer
     */
#ifdef MFEM_USE_MPI
    mfem::ParGridFunction* getModeShape(int modeIndex);
#else
    mfem::GridFunction* getModeShape(int modeIndex);
#endif

    /**
     * @brief Get the finite element space
     * @return Pointer to the Nedelec finite element space
     */
#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* getFiniteElementSpace() { return fespace_; }
#else
    mfem::FiniteElementSpace* getFiniteElementSpace() { return fespace_; }
#endif

private:
    /**
     * @brief Assembles the curl-curl stiffness matrix [A]
     * @details Implements the bilinear form: (∇×E, ∇×F)
     */
    void assembleCurlCurlMatrix();

    /**
     * @brief Assembles the mass matrix [M]
     * @details Implements the bilinear form: (E, F)
     */
    void assembleMassMatrix();

#ifdef MFEM_USE_MPI
    mfem::ParMesh* pmesh_;
    mfem::ParBilinearForm* curlCurlForm_;
    mfem::ParBilinearForm* massForm_;
    mfem::ParFiniteElementSpace* fespace_;
#else
    mfem::Mesh* mesh_;
    mfem::BilinearForm* curlCurlForm_;
    mfem::BilinearForm* massForm_;
    mfem::FiniteElementSpace* fespace_;
#endif

    // Nedelec (edge element) collection for vector fields
    mfem::FiniteElementCollection* feCollection_;
    
    // Stored eigenvalues and eigenvectors
    std::vector<double> eigenvalues_;
    std::vector<mfem::Vector> eigenvectors_;
    
    int polynomialOrder_;
    bool isSolved_;
};

} // namespace hpcfem

#endif // HPCFEM_PHYSICS_WAVEGUIDE_EIGEN_HPP
