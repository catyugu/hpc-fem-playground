/**
 * @file test_solver_hypre_amg.cpp
 * @brief Test for HypreAmgSolver using the Poisson MMS problem
 * 
 * Phase: 2, Step: 2.2
 * 
 * This test verifies that the HypreAmgSolver wrapper correctly solves
 * the Poisson problem using the SolverInterface abstraction.
 */

#include "../src/hpcfem/solver_hypre_amg.hpp"
#include "mfem.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;
using namespace hpcfem;

// Constants (reused from electrostatics test)
constexpr double PI = 3.141592653589793;
constexpr double MMS_TOLERANCE = 1.0e-3;
constexpr double CONDUCTIVITY = 1.0;
constexpr int MESH_ELEMENTS_1D = 8;
constexpr int POLYNOMIAL_ORDER = 2;

// Manufactured solution: u(x,y) = sin(πx)cos(πy)
double manufacturedSolution(const Vector &x)
{
    return std::sin(PI * x(0)) * std::cos(PI * x(1));
}

// Forcing term: f = 2π²sin(πx)cos(πy)
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
        std::cout << "Testing HypreAmgSolver with " << size << " DOFs" << std::endl;
    }

    // Setup the bilinear form
    ConstantCoefficient sigma(CONDUCTIVITY);
    
#ifdef MFEM_USE_MPI
    ParBilinearForm a(&fespace);
#else
    BilinearForm a(&fespace);
#endif
    a.AddDomainIntegrator(new DiffusionIntegrator(sigma));
    a.Assemble();

    // Setup the linear form
    FunctionCoefficient f_coeff(forcingTerm);
    
#ifdef MFEM_USE_MPI
    ParLinearForm b(&fespace);
#else
    LinearForm b(&fespace);
#endif
    b.AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
    b.Assemble();

    // Define solution vector with Dirichlet BC
    FunctionCoefficient u_exact_coeff(manufacturedSolution);
    
#ifdef MFEM_USE_MPI
    ParGridFunction x(&fespace);
#else
    GridFunction x(&fespace);
#endif
    x.ProjectCoefficient(u_exact_coeff);

    // Apply Dirichlet boundary conditions
    Array<int> ess_tdof_list;
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;
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

    // Create HypreAmgSolver using SolverInterface pointer
    HypreAmgSolver amg_solver(1e-12, 2000, 0);
    SolverInterface* solver = &amg_solver;

    // Solve using the interface
    solver->solve(A, B, X);

    // Recover the solution
    a.RecoverFEMSolution(X, b, x);

    // Compute L2 error norm
    double error = x.ComputeL2Error(u_exact_coeff);

    if (myid == 0)
    {
        std::cout << "L2 error: " << error << std::endl;
        std::cout << "Iterations: " << amg_solver.getNumIterations() << std::endl;
        std::cout << "Final norm: " << amg_solver.getFinalNorm() << std::endl;
    }

    // Test assertion
    bool test_passed = (error < MMS_TOLERANCE);

    // Cleanup
    delete pmesh;

    if (myid == 0)
    {
        if (test_passed)
        {
            std::cout << "Test passed: HypreAmgSolver solves Poisson problem" << std::endl;
        }
        else
        {
            std::cerr << "Test failed: L2 error " << error 
                      << " exceeds tolerance " << MMS_TOLERANCE << std::endl;
        }
    }

    return test_passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
