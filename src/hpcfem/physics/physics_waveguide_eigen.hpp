/**
 * @file physics/physics_waveguide_eigen.hpp
 * @brief 2D waveguide port eigenmode solver (TE/TM) for rectangular ports
 *
 * Implements Stage 1: 2D scalar Helmholtz eigenproblem for TE (Hz) and TM (Ez)
 * modes on a port cross-section. Supports serial (H1 + SLEPc) and parallel
 * (ParFiniteElementSpace + HypreLOBPCG) paths depending on MFEM_USE_MPI.
 */

#ifndef HPCFEM_PHYSICS_WAVEGUIDE_EIGEN_HPP
#define HPCFEM_PHYSICS_WAVEGUIDE_EIGEN_HPP

#include "mfem.hpp"
#include <vector>

namespace hpcfem
{

enum class WaveguideModeType { TE, TM };

/**
 * @brief 2D waveguide eigenmode solver for port cross-sections
 *
 * Solves the scalar Helmholtz eigenproblem on a 2D cross-section:
 *   -\nabla_t^2 psi = k_c^2 psi
 * with either Dirichlet (TM) or natural/Neumann (TE) boundary conditions.
 */
class PhysicsWaveguideEigen
{
public:
#ifdef MFEM_USE_MPI
    PhysicsWaveguideEigen(mfem::ParMesh *pmesh, int order, WaveguideModeType mode);
#else
    PhysicsWaveguideEigen(mfem::Mesh *mesh, int order, WaveguideModeType mode);
#endif

    ~PhysicsWaveguideEigen();

    /**
     * @brief Solve for the smallest `numModes` eigenvalues k_c^2
     * @param numModes Number of eigenvalues/modes to compute
     * @return vector of eigenvalues (k_c^2)
     */
    std::vector<double> solveEigenvalues(int numModes);

    /**
     * @brief Return local copies of eigenvectors (true/local DOF vectors)
     */
    std::vector<mfem::Vector> getEigenvectors() const;

    /**
     * @brief Get the computed mode as a GridFunction (caller deletes)
     */
#ifdef MFEM_USE_MPI
    mfem::ParGridFunction* getModeShape(int modeIndex);
#else
    mfem::GridFunction* getModeShape(int modeIndex);
#endif

private:
    void assembleStiffness();
    void assembleMass();

#ifdef MFEM_USE_MPI
    mfem::ParMesh* pmesh_ = nullptr;
    mfem::ParBilinearForm* stiffForm_ = nullptr;
    mfem::ParBilinearForm* massForm_ = nullptr;
    mfem::ParFiniteElementSpace* fespace_ = nullptr;
#else
    mfem::Mesh* mesh_ = nullptr;
    mfem::BilinearForm* stiffForm_ = nullptr;
    mfem::BilinearForm* massForm_ = nullptr;
    mfem::FiniteElementSpace* fespace_ = nullptr;
#endif

    mfem::FiniteElementCollection* fec_ = nullptr;
    WaveguideModeType modeType_;
    int order_ = 1;

    std::vector<double> eigenvalues_;
    std::vector<mfem::Vector> eigenvectors_;
    bool isSolved_ = false;
};

} // namespace hpcfem

#endif // HPCFEM_PHYSICS_WAVEGUIDE_EIGEN_HPP
