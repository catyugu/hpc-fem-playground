/**
 * @file physics/physics_joule_heating.hpp
 * @brief Coupled Joule heating physics implementation with nonlinear coupling
 */

#ifndef HPCFEM_PHYSICS_JOULE_HEATING_HPP
#define HPCFEM_PHYSICS_JOULE_HEATING_HPP

#include "hpcfem/core/physics_interface.hpp"
#include "mfem.hpp"
#include "joule_heating_coefficient.hpp"
#include <memory>

namespace hpcfem
{

class JouleHeatingPhysics
{
public:
#ifdef MFEM_USE_MPI
    JouleHeatingPhysics(mfem::ParMesh* pmesh,
                        int polynomialOrder,
                        double electricalConductivity,
                        double thermalConductivity);
#else
    JouleHeatingPhysics(mfem::Mesh* mesh,
                        int polynomialOrder,
                        double electricalConductivity,
                        double thermalConductivity);
#endif

    ~JouleHeatingPhysics();

#ifdef MFEM_USE_MPI
    void assemble(mfem::BlockOperator& blockOperator,
                 mfem::BlockVector& blockRHS,
                 mfem::BlockVector& blockSolution,
                 mfem::Array<int>& essTdofElectric,
                 mfem::Array<int>& essTdofThermal,
                 mfem::ParGridFunction* voltageGF = nullptr);

    mfem::ParFiniteElementSpace* getElectricSpace();
    mfem::ParFiniteElementSpace* getThermalSpace();
#else
    void assemble(mfem::BlockOperator& blockOperator,
                 mfem::BlockVector& blockRHS,
                 mfem::BlockVector& blockSolution,
                 mfem::Array<int>& essTdofElectric,
                 mfem::Array<int>& essTdofThermal,
                 mfem::GridFunction* voltageGF = nullptr);

    mfem::FiniteElementSpace* getElectricSpace();
    mfem::FiniteElementSpace* getThermalSpace();
#endif

    const mfem::Array<int>& getBlockOffsets() const;

private:
    mfem::H1_FECollection* fecElectric_;
    mfem::H1_FECollection* fecThermal_;
#ifdef MFEM_USE_MPI
    mfem::ParFiniteElementSpace* fespaceElectric_;
    mfem::ParFiniteElementSpace* fespaceThermal_;
    mfem::ParBilinearForm* bilinearElectric_;
    mfem::ParBilinearForm* bilinearThermal_;
    mfem::ParLinearForm* linearElectric_;
    mfem::ParLinearForm* linearThermal_;
    mfem::HypreParMatrix* matrixK_e_;
    mfem::HypreParMatrix* matrixK_t_;
    mfem::HypreParMatrix* matrixC_;
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

    double sigma_;
    double kappa_;
    mfem::ConstantCoefficient* sigmaCoeff_;
    mfem::ConstantCoefficient* kappaCoeff_;
    mfem::Array<int> blockOffsets_;
};

} // namespace hpcfem

#endif // HPCFEM_PHYSICS_JOULE_HEATING_HPP
