/**
 * @file test_problem_abstraction.cpp
 * @brief Test for the complete problem abstraction layer
 * 

 * 
 * This test validates the end-to-end workflow using:
 * - PhysicsInterface / ElectrostaticsPhysics
 * - SolverInterface / HypreAmgSolver
 * - FemProblem orchestrator
 */

#include "hpcfem/physics/physics_electrostatics.hpp"
#include "hpcfem/solvers/solver_hypre_amg.hpp"
#include "hpcfem/core/fem_problem.hpp"
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

#ifdef MFEM_USE_MPI
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    HYPRE_BigInt size = pmesh->GetGlobalNE();
#else
    Mesh *pmesh = mesh;
    int size = pmesh->GetNE();
#endif

    if (myid == 0)
    {
        std::cout << "Testing problem abstraction with " << size << " elements" << std::endl;
    }

    // Create coefficient objects
    ConstantCoefficient conductivity(CONDUCTIVITY);
    FunctionCoefficient forcing(forcingTerm);
    FunctionCoefficient boundary(manufacturedSolution);

    // Create physics module
    ElectrostaticsPhysics physics(pmesh, POLYNOMIAL_ORDER, &conductivity, &forcing, &boundary);

    // Create solver
    HypreAmgSolver solver(1e-12, 2000, 0);

    // Create problem orchestrator
    FemProblem problem(pmesh, &physics, &solver);

    // Assemble the system
    problem.assemble();

    if (myid == 0)
    {
        std::cout << "System assembled successfully" << std::endl;
    }

    // Solve the system
    problem.solve();

    if (myid == 0)
    {
        std::cout << "System solved successfully" << std::endl;
        std::cout << "Iterations: " << solver.getNumIterations() << std::endl;
        std::cout << "Final norm: " << solver.getFinalNorm() << std::endl;
    }

    // Compute L2 error norm
    auto* solutionGf = problem.getSolutionGridFunction();
    FunctionCoefficient exactSolution(manufacturedSolution);
    double error = solutionGf->ComputeL2Error(exactSolution);

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
            std::cout << "Test passed: Problem abstraction layer works correctly" << std::endl;
        }
        else
        {
            std::cerr << "Test failed: L2 error " << error 
                      << " exceeds tolerance " << MMS_TOLERANCE << std::endl;
        }
    }

    return test_passed ? EXIT_SUCCESS : EXIT_FAILURE;
}
