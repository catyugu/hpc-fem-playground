/**
 * @file test_ddm_one_level.cpp
 * @brief Test for one-level overlapping Schwarz (Additive Schwarz) preconditioner
 * 
 * Phase: 2, Step: 2.1
 * 
 * This test verifies correctness of the one-level domain decomposition method
 * on a 3D Poisson problem with manufactured solution.
 * 
 * Problem: -∇²u = f in Ω = [0,1]³
 *          u = u_exact on ∂Ω
 * 
 * Manufactured solution: u(x,y,z) = sin(πx) * sin(πy) * sin(πz)
 * Forcing term: f = 3π²sin(πx)sin(πy)sin(πz)
 */

#include "hpcfem/solver_one_level_schwarz.hpp"
#include "mfem.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;
using namespace hpcfem;

// Constants (Guideline #15: No magic numbers)
constexpr double PI = 3.141592653589793;
constexpr double MMS_TOLERANCE = 2.0e-3;  // L2 error tolerance (relaxed for coarse mesh)
constexpr double CONDUCTIVITY = 1.0;       // Diffusion coefficient
constexpr int MESH_ELEMENTS_1D = 4;        // Elements per dimension (small for testing)
constexpr int POLYNOMIAL_ORDER = 2;        // FE polynomial order
constexpr double CG_REL_TOL = 1.0e-12;
constexpr int CG_MAX_ITER = 2000;

// Manufactured solution: u(x,y,z) = sin(πx)sin(πy)sin(πz)
double manufacturedSolution(const Vector &x)
{
    return std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2));
}

// Forcing term: f = -∇²u = 3π²sin(πx)sin(πy)sin(πz)
double forcingTerm(const Vector &x)
{
    return 3.0 * PI * PI * std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2));
}

int main(int argc, char *argv[])
{
    // Initialize MPI (required for DDM)
    MPI_Session mpi(argc, argv);
    int myid = mpi.WorldRank();
    int numProcs = mpi.WorldSize();

    if (myid == 0)
    {
        std::cout << "Testing OneLevelSchwarz with " << numProcs << " MPI processes" << std::endl;
    }

    // Create a 3D Cartesian mesh [0,1]³
    Mesh *mesh = new Mesh(Mesh::MakeCartesian3D(
        MESH_ELEMENTS_1D, MESH_ELEMENTS_1D, MESH_ELEMENTS_1D, 
        Element::HEXAHEDRON));
    int dim = mesh->Dimension();

    // Partition mesh for parallel processing
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    // Define H1 finite element space
    H1_FECollection fec(POLYNOMIAL_ORDER, dim);
    ParFiniteElementSpace fespace(pmesh, &fec);
    HYPRE_BigInt size = fespace.GlobalTrueVSize();

    if (myid == 0)
    {
        std::cout << "Number of unknowns: " << size << std::endl;
    }

    // Setup the bilinear form a(u,v) = ∫ ∇u·∇v dx
    ConstantCoefficient sigma(CONDUCTIVITY);
    ParBilinearForm a(&fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(sigma));
    a.Assemble();

    // Setup the linear form b(v) = ∫ f*v dx
    FunctionCoefficient f_coeff(forcingTerm);
    ParLinearForm b(&fespace);
    b.AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
    b.Assemble();

    // Define solution vector with Dirichlet BC from manufactured solution
    FunctionCoefficient u_exact_coeff(manufacturedSolution);
    ParGridFunction x(&fespace);
    x.ProjectCoefficient(u_exact_coeff);

    // Apply Dirichlet boundary conditions (all boundaries)
    Array<int> ess_tdof_list;
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;  // All boundaries are essential
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // Form the linear system
    HypreParMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    // Create one-level Schwarz preconditioner
    OneLevelSchwarz prec(A);

    // Solve using CG with Schwarz preconditioner
    CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(CG_REL_TOL);
    cg.SetMaxIter(CG_MAX_ITER);
    cg.SetPrintLevel(myid == 0 ? 1 : 0);
    cg.SetPreconditioner(prec);
    cg.SetOperator(A);
    cg.Mult(B, X);

    if (myid == 0)
    {
        std::cout << "CG iterations: " << cg.GetNumIterations() << std::endl;
        std::cout << "CG final norm: " << cg.GetFinalNorm() << std::endl;
    }

    // Recover the solution
    a.RecoverFEMSolution(X, b, x);

    // Compute L2 error norm
    double error = x.ComputeL2Error(u_exact_coeff);

    if (myid == 0)
    {
        std::cout << "L2 error: " << error << std::endl;
    }

    // Test assertion: correctness check
    bool test_passed = (error < MMS_TOLERANCE) && (cg.GetConverged());

    // Cleanup
    delete pmesh;

    if (myid == 0)
    {
        if (test_passed)
        {
            std::cout << "Test passed: OneLevelSchwarz correctness test" << std::endl;
        }
        else
        {
            std::cerr << "Test failed: L2 error " << error 
                      << " exceeds tolerance " << MMS_TOLERANCE 
                      << " or CG did not converge" << std::endl;
        }
    }

    return test_passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
