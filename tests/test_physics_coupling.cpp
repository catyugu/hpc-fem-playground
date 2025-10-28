/**
 * @file test_physics_coupling.cpp
 * @brief Test for one-way coupled Joule heating simulation
 * 
 * Phase: 1, Step: 1.3
 * 
 * COUPLING:
 * 1. Solve electrostatics: -∇·(σ∇φ) = 0
 * 2. Compute heat source: Q = σ|∇φ|²
 * 3. Solve thermal: -∇·(k∇T) = Q
 * 
 * SANITY CHECKS:
 * - Maximum temperature > ambient temperature
 * - Maximum temperature < reasonable upper bound
 * - Heat generation is non-zero
 */

#include "mfem.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;

// Constants (Guideline #15: No magic numbers)
constexpr double ELECTRICAL_CONDUCTIVITY = 1.0e3;  // σ [S/m] - moderate conductor
constexpr double THERMAL_CONDUCTIVITY = 10.0;      // k [W/(m·K)] - better heat dissipation
constexpr double AMBIENT_TEMPERATURE = 293.15;     // T_amb [K]
constexpr double MAX_EXPECTED_TEMP = 400.0;        // Reasonable upper bound [K]
constexpr double VOLTAGE_HIGH = 1.0;               // High voltage boundary [V]
constexpr double VOLTAGE_LOW = 0.0;                // Low voltage boundary [V]
constexpr int MESH_ELEMENTS_1D = 10;
constexpr int POLYNOMIAL_ORDER = 2;
constexpr double MIN_HEAT_GENERATION = 1.0e-6;     // Minimum expected Q

// Joule heating coefficient: Q = σ|∇φ|²
class JouleHeatingCoefficient : public Coefficient
{
private:
    const GridFunction &potential_;
    double sigma_;
    mutable Vector grad_phi_;

public:
    JouleHeatingCoefficient(const GridFunction &phi, double sigma)
        : potential_(phi), sigma_(sigma) {}

    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        // Get gradient of potential: ∇φ
        potential_.GetGradient(T, grad_phi_);
        
        // Compute |∇φ|²
        double grad_phi_sq = grad_phi_ * grad_phi_;
        
        // Q = σ|∇φ|²
        return sigma_ * grad_phi_sq;
    }
};

int main(int argc, char *argv[])
{
    // Initialize MPI
#ifdef MFEM_USE_MPI
    MPI_Session mpi(argc, argv);
    int myid = mpi.WorldRank();
#else
    int myid = 0;
#endif

    // Create a 2D rectangular mesh [0,1] × [0,0.5]
    Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(
        MESH_ELEMENTS_1D, MESH_ELEMENTS_1D / 2, Element::QUADRILATERAL,
        false, 1.0, 0.5));
    int dim = mesh->Dimension();

#ifdef MFEM_USE_MPI
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
#else
    Mesh *pmesh = mesh;
#endif

    // ========================================================================
    // STEP 1: Solve Electrostatics Problem
    // ========================================================================
    
    if (myid == 0)
    {
        std::cout << "=== Step 1: Solving Electrostatics ===" << std::endl;
    }

#ifdef MFEM_USE_MPI
    H1_FECollection fec_elec(POLYNOMIAL_ORDER, dim);
    ParFiniteElementSpace fespace_elec(pmesh, &fec_elec);
#else
    H1_FECollection fec_elec(POLYNOMIAL_ORDER, dim);
    FiniteElementSpace fespace_elec(pmesh, &fec_elec);
#endif

    // Bilinear form: a(φ,v) = ∫ σ∇φ·∇v dx
    ConstantCoefficient sigma(ELECTRICAL_CONDUCTIVITY);
    
#ifdef MFEM_USE_MPI
    ParBilinearForm a_elec(&fespace_elec);
#else
    BilinearForm a_elec(&fespace_elec);
#endif
    a_elec.AddDomainIntegrator(new DiffusionIntegrator(sigma));
    a_elec.Assemble();

    // Linear form (zero RHS for Laplace equation)
#ifdef MFEM_USE_MPI
    ParLinearForm b_elec(&fespace_elec);
#else
    LinearForm b_elec(&fespace_elec);
#endif
    b_elec.Assemble();

    // Solution vector
#ifdef MFEM_USE_MPI
    ParGridFunction phi(&fespace_elec);
#else
    GridFunction phi(&fespace_elec);
#endif
    phi = 0.0;

    // Apply Dirichlet BCs: left = VOLTAGE_HIGH, right = VOLTAGE_LOW
    Array<int> ess_tdof_list_elec;
    Array<int> ess_bdr_elec(pmesh->bdr_attributes.Max());
    ess_bdr_elec = 0;
    ess_bdr_elec[0] = 1;  // Left boundary (attribute 1)
    ess_bdr_elec[1] = 1;  // Right boundary (attribute 2)
    fespace_elec.GetEssentialTrueDofs(ess_bdr_elec, ess_tdof_list_elec);

    // Set boundary values
    Array<int> ess_bdr_left(pmesh->bdr_attributes.Max());
    ess_bdr_left = 0;
    ess_bdr_left[0] = 1;  // Left boundary
    ConstantCoefficient voltage_high(VOLTAGE_HIGH);
    phi.ProjectBdrCoefficient(voltage_high, ess_bdr_left);

    Array<int> ess_bdr_right(pmesh->bdr_attributes.Max());
    ess_bdr_right = 0;
    ess_bdr_right[1] = 1;  // Right boundary
    ConstantCoefficient voltage_low(VOLTAGE_LOW);
    phi.ProjectBdrCoefficient(voltage_low, ess_bdr_right);

    // Form linear system
#ifdef MFEM_USE_MPI
    HypreParMatrix A_elec;
    Vector B_elec, X_elec;
    a_elec.FormLinearSystem(ess_tdof_list_elec, phi, b_elec, A_elec, X_elec, B_elec);
#else
    SparseMatrix A_elec;
    Vector B_elec, X_elec;
    a_elec.FormLinearSystem(ess_tdof_list_elec, phi, b_elec, A_elec, X_elec, B_elec);
#endif

    // Solve
#ifdef MFEM_USE_MPI
    HypreSolver *prec_elec = new HypreBoomerAMG(A_elec);
    CGSolver cg_elec(MPI_COMM_WORLD);
#else
    GSSmoother prec_elec(A_elec);
    CGSolver cg_elec;
#endif
    cg_elec.SetRelTol(1e-12);
    cg_elec.SetMaxIter(2000);
    cg_elec.SetPrintLevel(0);
    cg_elec.SetPreconditioner(*prec_elec);
    cg_elec.SetOperator(A_elec);
    cg_elec.Mult(B_elec, X_elec);

#ifdef MFEM_USE_MPI
    delete prec_elec;
#endif

    a_elec.RecoverFEMSolution(X_elec, b_elec, phi);

    if (myid == 0)
    {
        std::cout << "Electrostatics solved. Potential range: [" 
                  << phi.Min() << ", " << phi.Max() << "]" << std::endl;
    }

    // ========================================================================
    // STEP 2: Compute Joule Heating Source
    // ========================================================================
    
    if (myid == 0)
    {
        std::cout << "=== Step 2: Computing Joule Heating ===" << std::endl;
    }

    JouleHeatingCoefficient Q_coeff(phi, ELECTRICAL_CONDUCTIVITY);

    // Integrate to get total heat generation
#ifdef MFEM_USE_MPI
    ParLinearForm total_heat(&fespace_elec);
#else
    LinearForm total_heat(&fespace_elec);
#endif
    total_heat.AddDomainIntegrator(new DomainLFIntegrator(Q_coeff));
    total_heat.Assemble();
    double total_Q = total_heat.Sum();

    if (myid == 0)
    {
        std::cout << "Total heat generation: " << total_Q << " W" << std::endl;
    }

    // ========================================================================
    // STEP 3: Solve Thermal Problem
    // ========================================================================
    
    if (myid == 0)
    {
        std::cout << "=== Step 3: Solving Thermal Problem ===" << std::endl;
    }

#ifdef MFEM_USE_MPI
    H1_FECollection fec_therm(POLYNOMIAL_ORDER, dim);
    ParFiniteElementSpace fespace_therm(pmesh, &fec_therm);
#else
    H1_FECollection fec_therm(POLYNOMIAL_ORDER, dim);
    FiniteElementSpace fespace_therm(pmesh, &fec_therm);
#endif

    // Bilinear form: a(T,v) = ∫ k∇T·∇v dx
    ConstantCoefficient k(THERMAL_CONDUCTIVITY);
    
#ifdef MFEM_USE_MPI
    ParBilinearForm a_therm(&fespace_therm);
#else
    BilinearForm a_therm(&fespace_therm);
#endif
    a_therm.AddDomainIntegrator(new DiffusionIntegrator(k));
    a_therm.Assemble();

    // Linear form: b(v) = ∫ Q*v dx
#ifdef MFEM_USE_MPI
    ParLinearForm b_therm(&fespace_therm);
#else
    LinearForm b_therm(&fespace_therm);
#endif
    b_therm.AddDomainIntegrator(new DomainLFIntegrator(Q_coeff));
    b_therm.Assemble();

    // Solution vector
#ifdef MFEM_USE_MPI
    ParGridFunction T(&fespace_therm);
#else
    GridFunction T(&fespace_therm);
#endif
    
    // Set boundary to ambient temperature
    ConstantCoefficient T_amb(AMBIENT_TEMPERATURE);
    T.ProjectCoefficient(T_amb);

    // Apply Dirichlet BCs (all boundaries at ambient)
    Array<int> ess_tdof_list_therm;
    Array<int> ess_bdr_therm(pmesh->bdr_attributes.Max());
    ess_bdr_therm = 1;  // All boundaries
    fespace_therm.GetEssentialTrueDofs(ess_bdr_therm, ess_tdof_list_therm);

    // Form linear system
#ifdef MFEM_USE_MPI
    HypreParMatrix A_therm;
    Vector B_therm, X_therm;
    a_therm.FormLinearSystem(ess_tdof_list_therm, T, b_therm, A_therm, X_therm, B_therm);
#else
    SparseMatrix A_therm;
    Vector B_therm, X_therm;
    a_therm.FormLinearSystem(ess_tdof_list_therm, T, b_therm, A_therm, X_therm, B_therm);
#endif

    // Solve
#ifdef MFEM_USE_MPI
    HypreSolver *prec_therm = new HypreBoomerAMG(A_therm);
    CGSolver cg_therm(MPI_COMM_WORLD);
#else
    GSSmoother prec_therm(A_therm);
    CGSolver cg_therm;
#endif
    cg_therm.SetRelTol(1e-12);
    cg_therm.SetMaxIter(2000);
    cg_therm.SetPrintLevel(0);
    cg_therm.SetPreconditioner(*prec_therm);
    cg_therm.SetOperator(A_therm);
    cg_therm.Mult(B_therm, X_therm);

#ifdef MFEM_USE_MPI
    delete prec_therm;
#endif

    a_therm.RecoverFEMSolution(X_therm, b_therm, T);

    double max_temp = T.Max();
    double min_temp = T.Min();

    if (myid == 0)
    {
        std::cout << "Temperature range: [" << min_temp << ", " << max_temp << "] K" << std::endl;
    }

    // ========================================================================
    // SANITY CHECKS
    // ========================================================================
    
    bool test_passed = true;

    // Check 1: Total heat generation is non-zero
    if (total_Q < MIN_HEAT_GENERATION)
    {
        if (myid == 0)
        {
            std::cerr << "Test failed: Heat generation " << total_Q 
                      << " is below minimum " << MIN_HEAT_GENERATION << std::endl;
        }
        test_passed = false;
    }

    // Check 2: Maximum temperature > ambient
    if (max_temp <= AMBIENT_TEMPERATURE)
    {
        if (myid == 0)
        {
            std::cerr << "Test failed: Max temperature " << max_temp 
                      << " is not greater than ambient " << AMBIENT_TEMPERATURE << std::endl;
        }
        test_passed = false;
    }

    // Check 3: Maximum temperature < reasonable upper bound
    if (max_temp > MAX_EXPECTED_TEMP)
    {
        if (myid == 0)
        {
            std::cerr << "Test failed: Max temperature " << max_temp 
                      << " exceeds upper bound " << MAX_EXPECTED_TEMP << std::endl;
        }
        test_passed = false;
    }

    // Cleanup
    delete pmesh;

    if (myid == 0)
    {
        if (test_passed)
        {
            std::cout << "Test passed: One-way Joule heating coupling" << std::endl;
        }
    }

    return test_passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
