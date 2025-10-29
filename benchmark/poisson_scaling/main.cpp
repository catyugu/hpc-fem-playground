/**
 * @file main.cpp
 * @brief Benchmark application for Poisson solver scaling studies
 * 
 * This benchmark measures performance of AMG solver:
 * - HypreAmgSolver (algebraic multigrid)
 * 
 * Usage:
 *   ./benchmark_poisson --solver amg --mesh-size 32 --order 2
 *   mpirun -np 4 ./benchmark_poisson --solver amg --mesh-size 64
 */

#include "hpcfem/fem_problem.hpp"
#include "hpcfem/physics/physics_electrostatics.hpp"
#include "hpcfem/solvers/solver_hypre_amg.hpp"
#include "mfem.hpp"
#include <iostream>
#include <cmath>
#include <cstring>
#include <chrono>

using namespace mfem;
using namespace hpcfem;

// Constants
constexpr double PI = 3.141592653589793;
constexpr double CONDUCTIVITY = 1.0;
constexpr int DEFAULT_MESH_SIZE = 16;
constexpr int DEFAULT_ORDER = 2;

// Manufactured solution for testing
double manufacturedSolution(const Vector &x)
{
    return std::sin(PI * x(0)) * std::cos(PI * x(1));
}

double forcingTerm(const Vector &x)
{
    return 2.0 * PI * PI * std::sin(PI * x(0)) * std::cos(PI * x(1));
}

// Simple command-line parser
struct BenchmarkConfig
{
    const char* solverType;
    int meshSize;
    int polynomialOrder;
    bool printHelp;
    
    BenchmarkConfig() 
        : solverType("amg"), 
          meshSize(DEFAULT_MESH_SIZE), 
          polynomialOrder(DEFAULT_ORDER),
          printHelp(false) {}
};

void printUsage(const char* progName)
{
    std::cout << "Usage: " << progName << " [OPTIONS]\n"
              << "\nOptions:\n"
              << "  --solver TYPE      Solver type: 'amg' (default: amg)\n"
              << "  --mesh-size N      Mesh elements per dimension (default: 16)\n"
              << "  --order N          Polynomial order (default: 2)\n"
              << "  --help             Print this help message\n"
              << "\nExample:\n"
              << "  " << progName << " --solver amg --mesh-size 32 --order 3\n"
              << "  mpirun -np 4 " << progName << " --solver amg --mesh-size 64\n";
}

bool parseArgs(int argc, char* argv[], BenchmarkConfig& config)
{
    for (int i = 1; i < argc; ++i)
    {
        if (std::strcmp(argv[i], "--solver") == 0 && i + 1 < argc)
        {
            config.solverType = argv[++i];
        }
        else if (std::strcmp(argv[i], "--mesh-size") == 0 && i + 1 < argc)
        {
            config.meshSize = std::atoi(argv[++i]);
        }
        else if (std::strcmp(argv[i], "--order") == 0 && i + 1 < argc)
        {
            config.polynomialOrder = std::atoi(argv[++i]);
        }
        else if (std::strcmp(argv[i], "--help") == 0)
        {
            config.printHelp = true;
            return true;
        }
        else
        {
            std::cerr << "Unknown option: " << argv[i] << std::endl;
            return false;
        }
    }
    
    // Validate solver type
    if (std::strcmp(config.solverType, "amg") != 0)
    {
        std::cerr << "Error: solver must be 'amg'" << std::endl;
        return false;
    }
    
    // Validate mesh size
    if (config.meshSize < 1)
    {
        std::cerr << "Error: mesh-size must be positive" << std::endl;
        return false;
    }
    
    // Validate polynomial order
    if (config.polynomialOrder < 1)
    {
        std::cerr << "Error: order must be positive" << std::endl;
        return false;
    }
    
    return true;
}

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
    MPI_Session mpi(argc, argv);
    int myid = mpi.WorldRank();
    int numProcs = mpi.WorldSize();
#else
    int myid = 0;
    int numProcs = 1;
#endif

    // Parse command-line arguments
    BenchmarkConfig config;
    if (!parseArgs(argc, argv, config))
    {
        if (myid == 0)
        {
            printUsage(argv[0]);
        }
        return EXIT_FAILURE;
    }
    
    if (config.printHelp)
    {
        if (myid == 0)
        {
            printUsage(argv[0]);
        }
        return EXIT_SUCCESS;
    }

    // Print configuration
    if (myid == 0)
    {
        std::cout << "=== Poisson Solver Benchmark ===" << std::endl;
        std::cout << "Solver type: " << config.solverType << std::endl;
        std::cout << "Mesh size: " << config.meshSize << "x" << config.meshSize << std::endl;
        std::cout << "Polynomial order: " << config.polynomialOrder << std::endl;
        std::cout << "MPI processes: " << numProcs << std::endl;
    }

    // Create mesh
    Mesh *mesh = new Mesh(Mesh::MakeCartesian2D(
        config.meshSize, config.meshSize, Element::QUADRILATERAL));

#ifdef MFEM_USE_MPI
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
    HYPRE_BigInt globalDofs = 0;
#else
    Mesh *pmesh = mesh;
    int globalDofs = 0;
#endif

    // Create coefficients
    ConstantCoefficient conductivity(CONDUCTIVITY);
    FunctionCoefficient forcing(forcingTerm);
    FunctionCoefficient boundary(manufacturedSolution);

    // Create physics
    ElectrostaticsPhysics physics(pmesh, config.polynomialOrder, 
                                 &conductivity, &forcing, &boundary);

    // Get DOF count
#ifdef MFEM_USE_MPI
    globalDofs = physics.getFiniteElementSpace()->GlobalTrueVSize();
#else
    globalDofs = physics.getFiniteElementSpace()->GetTrueVSize();
#endif

    if (myid == 0)
    {
        std::cout << "Global DOFs: " << globalDofs << std::endl;
    }

    // Create solver based on configuration
    SolverInterface* solver = nullptr;
    HypreAmgSolver* amgSolver = nullptr;
    
    amgSolver = new HypreAmgSolver(1e-12, 2000, 0);
    solver = amgSolver;
    if (myid == 0)
    {
        std::cout << "Using AMG solver" << std::endl;
    }

    // Create problem
    FemProblem problem(pmesh, &physics, solver);

    // Assemble (timed)
    auto assembleStart = std::chrono::high_resolution_clock::now();
    problem.assemble();
    auto assembleEnd = std::chrono::high_resolution_clock::now();
    double assembleTime = std::chrono::duration<double>(assembleEnd - assembleStart).count();

    // Solve (timed)
    auto solveStart = std::chrono::high_resolution_clock::now();
    problem.solve();
    auto solveEnd = std::chrono::high_resolution_clock::now();
    double solveTime = std::chrono::duration<double>(solveEnd - solveStart).count();

    // Get iteration count
    int iterations = amgSolver->getNumIterations();

    // Compute error
    auto* solutionGf = problem.getSolutionGridFunction();
    FunctionCoefficient exactSolution(manufacturedSolution);
    double l2Error = solutionGf->ComputeL2Error(exactSolution);

    // Output parseable results
    if (myid == 0)
    {
        std::cout << "\n=== Results ===" << std::endl;
        std::cout << "DOFS: " << globalDofs << std::endl;
        std::cout << "ASSEMBLE_TIME: " << assembleTime << std::endl;
        std::cout << "SOLVE_TIME: " << solveTime << std::endl;
        std::cout << "TOTAL_TIME: " << (assembleTime + solveTime) << std::endl;
        std::cout << "ITERATIONS: " << iterations << std::endl;
        std::cout << "L2_ERROR: " << l2Error << std::endl;
        std::cout << "PROCS: " << numProcs << std::endl;
    }

    // Cleanup
    delete solver;
    delete pmesh;

    return EXIT_SUCCESS;
}
