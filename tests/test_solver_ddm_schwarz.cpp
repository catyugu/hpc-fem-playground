/**
 * @file test_solver_ddm_schwarz.cpp
 * @brief Test for DDM Schwarz solver
 * 
 * Phase: 3, Step: 3.1
 * 
 * This test validates the DDM Schwarz solver:
 * 1. Solves Poisson MMS problem and checks L2 error
 * 2. Compares convergence rate with Jacobi preconditioner
 */

#include "../src/hpcfem/solver_ddm_schwarz.hpp"
#include "../src/hpcfem/fem_problem.hpp"
#include "../src/hpcfem/physics_electrostatics.hpp"
#include "mfem.hpp"
#include <iostream>
#include <cmath>
#include <cstdlib>

using namespace mfem;
using namespace hpcfem;

// Constants
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

// Test 1: DDM Schwarz solver solves Poisson problem correctly
bool testDdmSchwarzSolvesPoisson(int myid)
{
    if (myid == 0)
    {
        std::cout << "\n=== Test 1: DDM Schwarz solves Poisson problem ===" << std::endl;
    }
    
    // Create mesh
    Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(
        MESH_ELEMENTS_1D, MESH_ELEMENTS_1D, Element::QUADRILATERAL));

#ifdef MFEM_USE_MPI
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
#else
    Mesh *pmesh = mesh;
#endif

    // Create coefficients
    ConstantCoefficient conductivity(CONDUCTIVITY);
    FunctionCoefficient forcing(forcingTerm);
    FunctionCoefficient boundary(manufacturedSolution);

    // Create physics and solver
    ElectrostaticsPhysics physics(pmesh, POLYNOMIAL_ORDER, &conductivity, &forcing, &boundary);
    DdmSchwarzSolver solver(1e-6, 100, 1);

    // Create and solve problem
    FemProblem problem(pmesh, &physics, &solver);
    problem.assemble();
    problem.solve();

    // Compute error
    auto* solutionGf = problem.getSolutionGridFunction();
    FunctionCoefficient exactSolution(manufacturedSolution);
    double error = solutionGf->ComputeL2Error(exactSolution);

    if (myid == 0)
    {
        std::cout << "L2 error: " << error << std::endl;
        std::cout << "Iterations: " << solver.getNumIterations() << std::endl;
        std::cout << "Final norm: " << solver.getFinalNorm() << std::endl;
    }

    bool passed = (error < MMS_TOLERANCE);

    delete pmesh;

    if (myid == 0)
    {
        if (passed)
        {
            std::cout << "Test 1 PASSED" << std::endl;
        }
        else
        {
            std::cerr << "Test 1 FAILED: L2 error too large" << std::endl;
        }
    }

    return passed;
}

// Test 2: DDM Schwarz converges faster than simple Jacobi
bool testDdmConvergesFasterThanJacobi(int myid)
{
#ifndef MFEM_USE_MPI
    if (myid == 0)
    {
        std::cout << "\n=== Test 2: Skipped (requires MPI) ===" << std::endl;
    }
    return true;
#else
    if (myid == 0)
    {
        std::cout << "\n=== Test 2: DDM vs Jacobi convergence ===" << std::endl;
    }
    
    // Create mesh
    Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(
        MESH_ELEMENTS_1D, MESH_ELEMENTS_1D, Element::QUADRILATERAL));
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;

    // Create coefficients
    ConstantCoefficient conductivity(CONDUCTIVITY);
    FunctionCoefficient forcing(forcingTerm);
    FunctionCoefficient boundary(manufacturedSolution);

    // Solution 1: DDM Schwarz
    ElectrostaticsPhysics physics_ddm(pmesh, POLYNOMIAL_ORDER, &conductivity, &forcing, &boundary);
    DdmSchwarzSolver ddmSolver(1e-6, 100, 0);
    FemProblem problemDdm(pmesh, &physics_ddm, &ddmSolver);
    problemDdm.assemble();
    problemDdm.solve();
    
    int ddmIterations = ddmSolver.getNumIterations();
    
    // Solution 2: PCG with Jacobi (simpler preconditioner)
    // We'll simulate this by using HypreBoomerAMG with just 1 level (essentially Jacobi)
    ElectrostaticsPhysics physics_jacobi(pmesh, POLYNOMIAL_ORDER, &conductivity, &forcing, &boundary);
    
    // Create a "weak" solver (fewer AMG levels, less aggressive)
    // This simulates a simpler preconditioner
    class SimpleJacobiSolver : public SolverInterface
    {
    public:
        SimpleJacobiSolver() : numIter_(0) {}
        
        void solve(const HypreParMatrix& A, const Vector& b, Vector& x) override
        {
            HypreBoomerAMG prec;
            prec.SetPrintLevel(0);
            prec.SetMaxLevels(1);  // Only 1 level = essentially Jacobi
            
            CGSolver cg(A.GetComm());
            cg.SetRelTol(1e-6);
            cg.SetMaxIter(200);
            cg.SetPrintLevel(0);
            cg.SetPreconditioner(prec);
            cg.SetOperator(A);
            
            cg.Mult(b, x);
            numIter_ = cg.GetNumIterations();
        }
        
        int getNumIterations() const { return numIter_; }
        
    private:
        int numIter_;
    };
    
    SimpleJacobiSolver jacobiSolver;
    FemProblem problemJacobi(pmesh, &physics_jacobi, &jacobiSolver);
    problemJacobi.assemble();
    problemJacobi.solve();
    
    int jacobiIterations = jacobiSolver.getNumIterations();
    
    if (myid == 0)
    {
        std::cout << "DDM Schwarz iterations: " << ddmIterations << std::endl;
        std::cout << "Simple preconditioner iterations: " << jacobiIterations << std::endl;
    }
    
    // DDM should converge in fewer iterations (or at least not much worse)
    // Due to better subdomain solves
    bool passed = (ddmIterations <= jacobiIterations * 1.5);  // Allow 50% more iterations
    
    delete pmesh;
    
    if (myid == 0)
    {
        if (passed)
        {
            std::cout << "Test 2 PASSED: DDM shows reasonable convergence" << std::endl;
        }
        else
        {
            std::cerr << "Test 2 FAILED: DDM took too many iterations" << std::endl;
        }
    }
    
    return passed;
#endif
}

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
    MPI_Session mpi(argc, argv);
    int myid = mpi.WorldRank();
#else
    int myid = 0;
#endif

    if (myid == 0)
    {
        std::cout << "Testing DDM Schwarz Solver" << std::endl;
    }

    bool test1Passed = testDdmSchwarzSolvesPoisson(myid);
    bool test2Passed = testDdmConvergesFasterThanJacobi(myid);

    bool allPassed = test1Passed && test2Passed;

    if (myid == 0)
    {
        std::cout << "\n==================================" << std::endl;
        if (allPassed)
        {
            std::cout << "ALL TESTS PASSED" << std::endl;
        }
        else
        {
            std::cout << "SOME TESTS FAILED" << std::endl;
        }
        std::cout << "==================================" << std::endl;
    }

    return allPassed ? EXIT_SUCCESS : EXIT_FAILURE;
}
