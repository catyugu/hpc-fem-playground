/**
 * @file test_mfem_baseline.cpp
 * @brief Baseline test replicating MFEM Example 1 (Poisson equation)
 * 
 * This test validates the parallel build by solving a simple Poisson problem
 * using MFEM's built-in functionality. It serves as a "known good" baseline.
 * 
 * Problem: -Δu = 1 in Ω = [0,1]³
 *          u = 0 on ∂Ω
 * 
 * Phase: 0, Step: 0.2
 */

#include "mfem.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;

// Constants (Guideline #15: No magic numbers)
constexpr int MESH_REFINEMENT_LEVELS = 2;
constexpr int POLYNOMIAL_ORDER = 1;
constexpr double EXPECTED_ERROR_NORM = 1.0e-2;  // Tolerance for this coarse mesh

int main(int argc, char *argv[])
{
    // Initialize MPI
#ifdef MFEM_USE_MPI
    MPI_Session mpi(argc, argv);
    int num_procs = mpi.WorldSize();
    int myid = mpi.WorldRank();
#else
    int num_procs = 1;
    int myid = 0;
#endif

    // Create a simple 3D Cartesian mesh
    Mesh *mesh = new Mesh(Mesh::MakeCartesian3D(4, 4, 4, Element::HEXAHEDRON));
    int dim = mesh->Dimension();
    
    // Refine the mesh
    for (int lev = 0; lev < MESH_REFINEMENT_LEVELS; lev++)
    {
        mesh->UniformRefinement();
    }

#ifdef MFEM_USE_MPI
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
#else
    Mesh *pmesh = mesh;
#endif

    // Define finite element space
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

    // Setup the linear form b (RHS) with constant forcing f = 1
    ConstantCoefficient one(1.0);
    
#ifdef MFEM_USE_MPI
    ParLinearForm b(&fespace);
#else
    LinearForm b(&fespace);
#endif
    b.AddDomainIntegrator(new DomainLFIntegrator(one));
    b.Assemble();

    // Setup the bilinear form a (LHS) for -Δu
#ifdef MFEM_USE_MPI
    ParBilinearForm a(&fespace);
#else
    BilinearForm a(&fespace);
#endif
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    a.Assemble();

    // Define solution vector and set homogeneous Dirichlet BC
#ifdef MFEM_USE_MPI
    ParGridFunction x(&fespace);
#else
    GridFunction x(&fespace);
#endif
    x = 0.0;

    Array<int> ess_tdof_list;
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;  // All boundaries are essential (Dirichlet)
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

    // Solve the linear system using HYPRE's AMG preconditioner
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

    // Compute error (approximate check - we don't have exact solution for arbitrary mesh)
    // For this simple problem, we check that the solution has reasonable values
    double max_val = x.Max();
    double min_val = x.Min();
    
    if (myid == 0)
    {
        std::cout << "Solution range: [" << min_val << ", " << max_val << "]" << std::endl;
    }

    // Sanity checks
    constexpr double MIN_EXPECTED = 0.0;
    constexpr double MAX_EXPECTED = 1.0;
    bool test_passed = true;

    if (min_val < MIN_EXPECTED - EXPECTED_ERROR_NORM)
    {
        if (myid == 0)
        {
            std::cerr << "Test failed: Minimum value " << min_val 
                      << " is below expected range" << std::endl;
        }
        test_passed = false;
    }

    if (max_val > MAX_EXPECTED + EXPECTED_ERROR_NORM)
    {
        if (myid == 0)
        {
            std::cerr << "Test failed: Maximum value " << max_val 
                      << " exceeds expected range" << std::endl;
        }
        test_passed = false;
    }

    // Cleanup
    delete pmesh;

    if (myid == 0)
    {
        if (test_passed)
        {
            std::cout << "Test passed: MFEM baseline solution computed successfully" << std::endl;
        }
    }

    return test_passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
