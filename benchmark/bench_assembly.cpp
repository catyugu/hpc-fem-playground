#include "mfem.hpp"
#include "plugins/mfem_plugins.hpp"
#include <iostream>
#include <iomanip>

using namespace mfem;
using namespace mfem_plugins;

// Benchmark assembly strategies for Poisson problem
void benchmark_assembly(const char* mesh_file, int order, int refine_levels) {
    std::cout << "\n========================================" << std::endl;
    std::cout << "Assembly Strategy Benchmark" << std::endl;
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
    
    // Test vectors
    Vector x(fes.GetTrueVSize());
    Vector y(fes.GetTrueVSize());
    x.Randomize();
    
    // Coefficient for diffusion
    ConstantCoefficient one(1.0);
    
    // Results table header
    std::cout << std::setw(20) << "Strategy" << " | "
              << std::setw(12) << "Setup (ms)" << " | "
              << std::setw(12) << "Apply (ms)" << " | "
              << std::setw(12) << "Memory (MB)" << " | "
              << std::setw(10) << "Speedup" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    
    double baseline_time = 0.0;
    
    // Test each assembly strategy
    AssemblyStrategy strategies[] = {
        AssemblyStrategy::FULL,
        // Note: PARTIAL and MATRIX_FREE require special integrator support in MFEM
        // They are not universally available for all integrators
        // AssemblyStrategy::PARTIAL,
        // AssemblyStrategy::MATRIX_FREE
    };
    
    const char* strategy_names[] = {
        "Full Assembly",
        // "Partial Assembly",
        // "Matrix-Free"
    };
    
    for (size_t i = 0; i < sizeof(strategies) / sizeof(strategies[0]); i++) {
        // Create bilinear form
        BilinearForm a(&fes);
        a.AddDomainIntegrator(new DiffusionIntegrator(one));
        
        // Create plugin
        auto plugin = create_assembly_plugin(strategies[i]);
        plugin->set_fespace(&fes);
        
        // Assemble
        plugin->assemble(a);
        
        // Test matrix-vector product (multiple iterations for timing)
        Timer apply_timer;
        apply_timer.start();
        int num_applies = 10;
        for (int j = 0; j < num_applies; j++) {
            plugin->mult(x, y);
        }
        double apply_time = apply_timer.stop() / num_applies;
        
        // Get metrics
        double setup_time = plugin->metrics().setup_time_ms;
        double memory_mb = plugin->metrics().memory_bytes / (1024.0 * 1024.0);
        
        if (i == 0) baseline_time = setup_time + apply_time;
        double speedup = baseline_time / (setup_time + apply_time);
        
        std::cout << std::setw(20) << strategy_names[i] << " | "
                  << std::setw(12) << std::fixed << std::setprecision(3) << setup_time << " | "
                  << std::setw(12) << std::fixed << std::setprecision(3) << apply_time << " | "
                  << std::setw(12) << std::fixed << std::setprecision(2) << memory_mb << " | "
                  << std::setw(10) << std::fixed << std::setprecision(2) << speedup << "x" << std::endl;
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
    
    benchmark_assembly(mesh_file, order, refine_levels);
    
    return 0;
}
