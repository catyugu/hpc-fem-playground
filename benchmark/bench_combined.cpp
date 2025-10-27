#include "mfem.hpp"
#include "plugins/mfem_plugins.hpp"
#include <iostream>
#include <iomanip>

using namespace mfem;
using namespace mfem_plugins;

// Benchmark combined assembly + solver workflow
void benchmark_combined_workflow(const char* mesh_file, int order, int refine_levels) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Combined Workflow Benchmark" << std::endl;
    std::cout << "Assembly + Solver Combinations" << std::endl;
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
    
    // Setup problem coefficients
    ConstantCoefficient one(1.0);
    
    // Essential BCs
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1;
    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    
    // RHS
    LinearForm b(&fes);
    b = 0.0;
    b.Assemble();
    
    GridFunction x(&fes);
    x = 0.0;
    
    // Results table header
    std::cout << std::setw(20) << "Assembly" << " | "
              << std::setw(20) << "Solver" << " | "
              << std::setw(12) << "Asm (ms)" << " | "
              << std::setw(12) << "Solve (ms)" << " | "
              << std::setw(12) << "Total (ms)" << " | "
              << std::setw(8) << "Iters" << std::endl;
    std::cout << std::string(110, '-') << std::endl;
    
    // Test configurations
    struct Config {
        AssemblyStrategy asm_strategy;
        SolverType solver_type;
        PreconditionerType precond;
        const char* asm_name;
        const char* solver_name;
    };
    
    std::vector<Config> configs = {
        // Standard full assembly with various solvers
        {AssemblyStrategy::FULL, SolverType::CG, PreconditionerType::NONE, 
         "Full Assembly", "CG (no precond)"},
        {AssemblyStrategy::FULL, SolverType::CG, PreconditionerType::JACOBI, 
         "Full Assembly", "CG + Jacobi"},
        {AssemblyStrategy::FULL, SolverType::CG, PreconditionerType::GAUSS_SEIDEL, 
         "Full Assembly", "CG + GS"},
         
        // Note: PARTIAL, MATRIX_FREE, ELEMENT_BY_ELEMENT require special support
        // Commented out to avoid runtime errors
        /*
        // Partial assembly (element matrices)
        {AssemblyStrategy::PARTIAL, SolverType::CG, PreconditionerType::JACOBI, 
         "Partial Assembly", "CG + Jacobi"},
         
        // Matrix-free with preconditioned CG
        {AssemblyStrategy::MATRIX_FREE, SolverType::CG, PreconditionerType::JACOBI, 
         "Matrix-Free", "CG + Jacobi"},
         
        // Element-by-element
        {AssemblyStrategy::ELEMENT_BY_ELEMENT, SolverType::CG, PreconditionerType::JACOBI, 
         "Element-by-Element", "CG + Jacobi"},
        */
    };
    
    for (const auto& config : configs) {
        // Reset solution
        x = 0.0;
        
        // Create assembly plugin
        auto assembly = create_assembly_plugin(config.asm_strategy);
        
        // Assemble bilinear form
        BilinearForm a(&fes);
        a.AddDomainIntegrator(new DiffusionIntegrator(one));
        
        assembly->assemble(a);
        
        // Form system
        SparseMatrix A;
        Vector B, X;
        a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
        
        // Add some RHS for testing
        B.Randomize();
        
        // Create solver plugin
        auto solver = create_solver_plugin(config.solver_type, config.precond);
        solver->set_tolerance(1e-9, 1e-12);
        solver->set_max_iterations(1000);
        solver->set_print_level(0);
        
        // Setup and solve
        solver->setup(A);
        solver->solve(B, X);
        
        // Get metrics
        const auto& asm_metrics = assembly->metrics();
        const auto& solver_metrics = solver->metrics();
        
        double total_time = asm_metrics.setup_time_ms + solver_metrics.total_time_ms;
        
        std::cout << std::setw(20) << config.asm_name << " | "
                  << std::setw(20) << config.solver_name << " | "
                  << std::setw(12) << std::fixed << std::setprecision(3) << asm_metrics.setup_time_ms << " | "
                  << std::setw(12) << std::fixed << std::setprecision(3) << solver_metrics.solve_time_ms << " | "
                  << std::setw(12) << std::fixed << std::setprecision(3) << total_time << " | "
                  << std::setw(8) << solver_metrics.iterations << std::endl;
    }
    
    std::cout << "========================================\n" << std::endl;
    
    // Performance recommendations
    std::cout << "\nPerformance Recommendations:" << std::endl;
    std::cout << "- For small/medium problems: Full assembly + preconditioned CG" << std::endl;
    std::cout << "- For large problems with many matvecs: Matrix-free + Jacobi" << std::endl;
    std::cout << "- For memory-constrained systems: Partial/Element-by-element assembly" << std::endl;
    std::cout << "- For complex geometries: Full assembly + direct solver (if available)" << std::endl;
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
    
    benchmark_combined_workflow(mesh_file, order, refine_levels);
    
    return 0;
}
