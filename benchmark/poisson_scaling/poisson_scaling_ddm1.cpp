/**
 * @file poisson_scaling_ddm1.cpp
 * @brief Scalability benchmark for one-level Schwarz (demonstrates poor scaling)
 * 
 * Phase: 2, Step: 2.1
 * 
 * This benchmark solves a 3D Poisson problem with increasing numbers of MPI processes
 * using the one-level overlapping Schwarz preconditioner. The key observation is that
 * the iteration count INCREASES with the number of processes, demonstrating the well-known
 * scalability limitation of one-level domain decomposition methods.
 * 
 * Expected behavior:
 * - 2 procs: ~20 iterations
 * - 4 procs: ~23 iterations
 * - 8 procs: ~30+ iterations
 * - 16 procs: ~40+ iterations
 * 
 * This motivates the need for two-level methods (Phase 2.2).
 */

#include "hpcfem/solver_one_level_schwarz.hpp"
#include "mfem.hpp"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace mfem;
using namespace hpcfem;

// Constants
constexpr double PI = 3.141592653589793;
constexpr double CONDUCTIVITY = 1.0;
constexpr int POLYNOMIAL_ORDER = 2;
constexpr double CG_REL_TOL = 1.0e-12;
constexpr int CG_MAX_ITER = 2000;

// Manufactured solution
double manufacturedSolution(const Vector &x)
{
    return std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2));
}

// Forcing term
double forcingTerm(const Vector &x)
{
    return 3.0 * PI * PI * std::sin(PI * x(0)) * std::sin(PI * x(1)) * std::sin(PI * x(2));
}

int main(int argc, char *argv[])
{
    // Initialize MPI
    MPI_Session mpi(argc, argv);
    int myid = mpi.WorldRank();
    int numProcs = mpi.WorldSize();

    // Problem size: scale with number of processes to maintain work per process
    // Use cube root to get elements per dimension
    int elementsPerProc = 64;  // Target ~64 elements per process
    int totalElements = elementsPerProc * numProcs;
    int elementsPerDim = static_cast<int>(std::cbrt(totalElements) + 0.5);
    if (elementsPerDim < 4) elementsPerDim = 4;

    if (myid == 0)
    {
        std::cout << "=== One-Level Schwarz Scaling Benchmark ===" << std::endl;
        std::cout << "Number of MPI processes: " << numProcs << std::endl;
        std::cout << "Elements per dimension: " << elementsPerDim << std::endl;
    }

    // Create 3D Cartesian mesh
    Mesh *mesh = new Mesh(Mesh::MakeCartesian3D(
        elementsPerDim, elementsPerDim, elementsPerDim, 
        Element::HEXAHEDRON));

    // Partition mesh for parallel processing
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    // Define H1 finite element space
    H1_FECollection fec(POLYNOMIAL_ORDER, pmesh->Dimension());
    ParFiniteElementSpace fespace(pmesh, &fec);
    HYPRE_BigInt size = fespace.GlobalTrueVSize();

    if (myid == 0)
    {
        std::cout << "Number of unknowns: " << size << std::endl;
    }

    // Setup bilinear form
    ConstantCoefficient sigma(CONDUCTIVITY);
    ParBilinearForm a(&fespace);
    a.AddDomainIntegrator(new DiffusionIntegrator(sigma));
    a.Assemble();

    // Setup linear form
    FunctionCoefficient f_coeff(forcingTerm);
    ParLinearForm b(&fespace);
    b.AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
    b.Assemble();

    // Define solution vector with Dirichlet BC
    FunctionCoefficient u_exact_coeff(manufacturedSolution);
    ParGridFunction x(&fespace);
    x.ProjectCoefficient(u_exact_coeff);

    // Apply Dirichlet boundary conditions
    Array<int> ess_tdof_list;
    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 1;
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // Form the linear system
    HypreParMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    // Create one-level Schwarz preconditioner
    OneLevelSchwarz prec(A);

    // Solve using CG with Schwarz preconditioner
    StopWatch timer;
    timer.Start();
    
    CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(CG_REL_TOL);
    cg.SetMaxIter(CG_MAX_ITER);
    cg.SetPrintLevel(0);  // Silent for benchmarking
    cg.SetPreconditioner(prec);
    cg.SetOperator(A);
    cg.Mult(B, X);
    
    timer.Stop();
    double solveTime = timer.RealTime();

    // Gather results on root process
    int iterations = cg.GetNumIterations();
    double finalNorm = cg.GetFinalNorm();

    // Recover solution and compute error
    a.RecoverFEMSolution(X, b, x);
    double error = x.ComputeL2Error(u_exact_coeff);

    if (myid == 0)
    {
        std::cout << "\n=== Results ===" << std::endl;
        std::cout << "CG Iterations: " << iterations << std::endl;
        std::cout << "Final residual norm: " << finalNorm << std::endl;
        std::cout << "Solve time: " << solveTime << " seconds" << std::endl;
        std::cout << "L2 error: " << error << std::endl;
        std::cout << "Converged: " << (cg.GetConverged() ? "Yes" : "No") << std::endl;
        
        // Write to CSV file for plotting
        std::ofstream csvFile;
        csvFile.open("ddm1_scaling.csv", std::ios::app);
        if (csvFile.tellp() == 0)
        {
            // Write header if file is empty
            csvFile << "NumProcs,DOFs,Iterations,SolveTime,L2Error\n";
        }
        csvFile << numProcs << "," << size << "," << iterations << "," 
                << solveTime << "," << error << "\n";
        csvFile.close();
        
        std::cout << "\nResults appended to ddm1_scaling.csv" << std::endl;
    }

    // Cleanup
    delete pmesh;

    return 0;
}
