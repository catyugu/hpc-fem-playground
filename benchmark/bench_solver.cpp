#include "mfem.hpp"
#include "plugins/mfem_plugins.hpp"
#include <iostream>
#include <iomanip>

using namespace mfem;
using namespace mfem_plugins;

// Benchmark solver strategies for Poisson problem
void benchmark_solvers(const char* mesh_file, int order, int refine_levels) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Solver Strategy Benchmark" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Mesh: " << mesh_file << std::endl;
    std::cout << "Order: " << order << std::endl;
    std::cout << "Refinement levels: " << refine_levels << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Load mesh
    Mesh mesh(mesh_file);
    int dim = mesh.Dimension();
    
    // Refine mesh
    for (int i = 0; i < refine_levels; i++) {
        mesh.UniformRefinement();
    }
    
    std::cout << "Mesh elements: " << mesh.GetNE() << std::endl;
    
    // Create finite element space
    H1_FECollection fec(order, dim);
    FiniteElementSpace fes(&mesh, &fec);
    
    std::cout << "Number of DOFs: " << fes.GetTrueVSize() << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Build system: Poisson equation
    BilinearForm a(&fes);
    ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    
    LinearForm b(&fes);
    b = 0.0;
    b.Assemble();
    
    // Essential BCs (homogeneous Dirichlet)
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1;
    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    
    GridFunction x(&fes);
    x = 0.0;
    
    // Assemble and form system
    a.Assemble();
    a.Finalize();
    
    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
    
    // Add some non-zero RHS for testing
    B.Randomize();
    
    std::cout << "System size: " << A.Height() << " x " << A.Width() << std::endl;
    std::cout << "NNZ: " << A.NumNonZeroElems() << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Results table header
    std::cout << std::setw(25) << "Solver" << " | "
              << std::setw(12) << "Setup (ms)" << " | "
              << std::setw(12) << "Solve (ms)" << " | "
              << std::setw(10) << "Total (ms)" << " | "
              << std::setw(8) << "Iters" << " | "
              << std::setw(10) << "Residual" << std::endl;
    std::cout << std::string(100, '-') << std::endl;
    
    // Test configurations
    struct SolverConfig {
        SolverType type;
        PreconditionerType precond;
        const char* name;
    };
    
    std::vector<SolverConfig> configs = {
        {SolverType::CG, PreconditionerType::NONE, "CG (no precond)"},
        {SolverType::CG, PreconditionerType::JACOBI, "CG + Jacobi"},
        {SolverType::CG, PreconditionerType::GAUSS_SEIDEL, "CG + Gauss-Seidel"},
        {SolverType::GMRES, PreconditionerType::NONE, "GMRES (no precond)"},
        {SolverType::GMRES, PreconditionerType::JACOBI, "GMRES + Jacobi"},
    };
    
    // Add direct solver if available
#ifdef MFEM_USE_SUITESPARSE
    configs.push_back({SolverType::DIRECT, PreconditionerType::NONE, "Direct (UMFPACK)"});
#endif
    
    for (const auto& config : configs) {
        // Create plugin
        auto solver = create_solver_plugin(config.type, config.precond);
        solver->set_tolerance(1e-9, 1e-12);
        solver->set_max_iterations(1000);
        solver->set_print_level(0);
        
        // Setup solver
        solver->setup(A);
        
        // Solve
        Vector x_solve = X;  // Copy initial guess
        solver->solve(B, x_solve);
        
        // Get metrics
        const auto& metrics = solver->metrics();
        
        std::cout << std::setw(25) << config.name << " | "
                  << std::setw(12) << std::fixed << std::setprecision(3) << metrics.setup_time_ms << " | "
                  << std::setw(12) << std::fixed << std::setprecision(3) << metrics.solve_time_ms << " | "
                  << std::setw(10) << std::fixed << std::setprecision(3) << metrics.total_time_ms << " | "
                  << std::setw(8) << metrics.iterations << " | "
                  << std::setw(10) << std::scientific << std::setprecision(2) << metrics.residual << std::endl;
    }
    
    std::cout << "========================================\n" << std::endl;
}

int main(int argc, char *argv[]) {
    const char* mesh_file = "testdata/testmesh_cube.mesh";
    int order = 1;
    int refine_levels = 0;
    
    if (argc > 1) mesh_file = argv[1];
    if (argc > 2) order = atoi(argv[2]);
    if (argc > 3) refine_levels = atoi(argv[3]);
    
    print_plugin_info();
    
    benchmark_solvers(mesh_file, order, refine_levels);
    
    return 0;
}
