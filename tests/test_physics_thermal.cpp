/**
 * @file test_physics_thermal.cpp
 * @brief Test for standalone thermal solver using Method of Manufactured Solutions
 * 
 * Phase: 1, Step: 1.2
 * 
 * Problem: -∇·(k∇T) = Q in Ω = [0,1]²
 *          T = T_exact on ∂Ω
 * 
 * Manufactured solution: T(x,y) = cosh(x) * sinh(y)
 * Forcing term: Q = -∇²T = -k[sinh(x)sinh(y) + cosh(x)cosh(y)]
 */

#include "mfem.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;

// Constants (Guideline #15: No magic numbers)
constexpr double MMS_TOLERANCE = 1.0e-3;  // L2 error tolerance
constexpr double THERMAL_CONDUCTIVITY = 1.0;  // k = 1.0 [W/(m·K)]
constexpr int MESH_ELEMENTS_1D = 8;
constexpr int POLYNOMIAL_ORDER = 2;

// Manufactured solution: T(x,y) = cosh(x) * sinh(y)
double manufacturedTemperature(const Vector &x)
{
    return std::cosh(x(0)) * std::sinh(x(1));
}

// Forcing term: Q = -∇²T = -[∂²T/∂x² + ∂²T/∂y²]
// ∂²T/∂x² = cosh(x)sinh(y), ∂²T/∂y² = cosh(x)sinh(y)
// Q = -[cosh(x)sinh(y) + cosh(x)sinh(y)] = -2cosh(x)sinh(y)
double heatSource(const Vector &x)
{
    return -2.0 * std::cosh(x(0)) * std::sinh(x(1));
}

int main(int argc, char *argv[])
{
    // Initialize MPI
#ifdef MFEM_USE_MPI
    MPI_Session mpi(argc, argv);
    int myid = mpi.WorldRank();
#else
    int myid = 0;
#endif

    // Create a 2D Cartesian mesh [0,1]²
    Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(
        MESH_ELEMENTS_1D, MESH_ELEMENTS_1D, Element::QUADRILATERAL));
    int dim = mesh->Dimension();

#ifdef MFEM_USE_MPI
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
#else
    Mesh *pmesh = mesh;
#endif

    // Define H1 finite element space
#ifdef MFEM_USE_MPI
    H1_FECollection fec(POLYNOMIAL_ORDER, dim);
    ParFiniteElementSpace fespace(pmesh, &fec);
    HYPRE_BigInt size = fespace.GlobalTrueVSize();
#else
    H1_FECollection fec(POLYNOMIAL_ORDER, dim);
    FiniteElementSpace fespace(pmesh, &fec);
    int size = fespace.GetTrueVSize();
#endif

    if (myid == 0)
    {
        std::cout << "Number of unknowns: " << size << std::endl;
    }

    // Setup the bilinear form a(T,v) = ∫ k∇T·∇v dx
    ConstantCoefficient k(THERMAL_CONDUCTIVITY);
    
#ifdef MFEM_USE_MPI
    ParBilinearForm a(&fespace);
#else
    BilinearForm a(&fespace);
#endif
    a.AddDomainIntegrator(new DiffusionIntegrator(k));
    a.Assemble();

    // Setup the linear form b(v) = ∫ Q*v dx
    FunctionCoefficient Q_coeff(heatSource);
    
#ifdef MFEM_USE_MPI
    ParLinearForm b(&fespace);
#else
    LinearForm b(&fespace);
#endif
    b.AddDomainIntegrator(new DomainLFIntegrator(Q_coeff));
    b.Assemble();

    // Define solution vector with Dirichlet BC from manufactured solution
    FunctionCoefficient T_exact_coeff(manufacturedTemperature);
    
#ifdef MFEM_USE_MPI
    ParGridFunction T(&fespace);
#else
    GridFunction T(&fespace);
#endif
    T.ProjectCoefficient(T_exact_coeff);

    // Apply Dirichlet boundary conditions (all boundaries)
    Array<int> ess_tdof_list;
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;  // All boundaries are essential
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // Form the linear system
#ifdef MFEM_USE_MPI
    HypreParMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, T, b, A, X, B);
#else
    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, T, b, A, X, B);
#endif

    // Solve using HYPRE AMG preconditioner
#ifdef MFEM_USE_MPI
    HypreSolver *prec = new HypreBoomerAMG(A);
    CGSolver cg(MPI_COMM_WORLD);
#else
    GSSmoother prec(A);
    CGSolver cg;
#endif
    cg.SetRelTol(1e-12);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(0);
    cg.SetPreconditioner(*prec);
    cg.SetOperator(A);
    cg.Mult(B, X);

#ifdef MFEM_USE_MPI
    delete prec;
#endif

    // Recover the solution
    a.RecoverFEMSolution(X, b, T);

    // Compute L2 error norm
    double error = T.ComputeL2Error(T_exact_coeff);

    if (myid == 0)
    {
        std::cout << "L2 error: " << error << std::endl;
    }

    // Test assertion
    bool test_passed = (error < MMS_TOLERANCE);

    // Cleanup
    delete pmesh;

    if (myid == 0)
    {
        if (test_passed)
        {
            std::cout << "Test passed: Thermal MMS test" << std::endl;
        }
        else
        {
            std::cerr << "Test failed: L2 error " << error 
                      << " exceeds tolerance " << MMS_TOLERANCE << std::endl;
        }
    }

    return test_passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
