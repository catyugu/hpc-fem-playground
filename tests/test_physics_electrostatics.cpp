/**
 * @file test_physics_electrostatics.cpp
 * @brief Test for standalone electrostatics solver using Method of Manufactured Solutions
 * 

 * 
 * Problem: -∇·(σ∇u) = f in Ω = [0,1]²
 *          u = u_exact on ∂Ω
 * 
 * Manufactured solution: u(x,y) = sin(πx) * cos(πy)
 * Forcing term: f = -∇²u = 2π²sin(πx)cos(πy)
 */

#include "mfem.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;

// Constants (Guideline #15: No magic numbers)
constexpr double PI = 3.141592653589793;
constexpr double MMS_TOLERANCE = 1.0e-3;  // L2 error tolerance
constexpr double CONDUCTIVITY = 1.0;       // σ = 1.0
constexpr int MESH_ELEMENTS_1D = 8;
constexpr int POLYNOMIAL_ORDER = 2;

// Manufactured solution: u(x,y) = sin(πx)cos(πy)
double manufactureSolution(const Vector &x)
{
    return std::sin(PI * x(0)) * std::cos(PI * x(1));
}

// Forcing term: f = -∇²u = 2π²sin(πx)cos(πy)
double forcingTerm(const Vector &x)
{
    return 2.0 * PI * PI * std::sin(PI * x(0)) * std::cos(PI * x(1));
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

    // Setup the bilinear form a(u,v) = ∫ σ∇u·∇v dx
    ConstantCoefficient sigma(CONDUCTIVITY);
    
#ifdef MFEM_USE_MPI
    ParBilinearForm a(&fespace);
#else
    BilinearForm a(&fespace);
#endif
    a.AddDomainIntegrator(new DiffusionIntegrator(sigma));
    a.Assemble();

    // Setup the linear form b(v) = ∫ f*v dx
    FunctionCoefficient f_coeff(forcingTerm);
    
#ifdef MFEM_USE_MPI
    ParLinearForm b(&fespace);
#else
    LinearForm b(&fespace);
#endif
    b.AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
    b.Assemble();

    // Define solution vector with Dirichlet BC from manufactured solution
    FunctionCoefficient u_exact_coeff(manufactureSolution);
    
#ifdef MFEM_USE_MPI
    ParGridFunction x(&fespace);
#else
    GridFunction x(&fespace);
#endif
    x.ProjectCoefficient(u_exact_coeff);

    // Apply Dirichlet boundary conditions (all boundaries)
    Array<int> ess_tdof_list;
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;  // All boundaries are essential
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // Form the linear system
#ifdef MFEM_USE_MPI
    HypreParMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
#else
    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
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
    a.RecoverFEMSolution(X, b, x);

    // Compute L2 error norm
    double error = x.ComputeL2Error(u_exact_coeff);

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
            std::cout << "Test passed: Electrostatics MMS test" << std::endl;
        }
        else
        {
            std::cerr << "Test failed: L2 error " << error 
                      << " exceeds tolerance " << MMS_TOLERANCE << std::endl;
        }
    }

    return test_passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
