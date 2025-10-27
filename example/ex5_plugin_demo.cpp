#include "mfem.hpp"
#include "plugins/mfem_plugins.hpp"
#include <iostream>
#include <iomanip>

using namespace mfem;
using namespace mfem_plugins;

// Example demonstrating plugin-based workflow customization
// This shows how to easily swap assembly and solver strategies

int main(int argc, char *argv[]) {
    // Parse command line arguments
    const char* mesh_file = "testdata/testmesh_cube.mesh";
    int order = 2;
    int refine_levels = 1;
    
    // Plugin configuration (can be changed at runtime)
    AssemblyStrategy asm_strategy = AssemblyStrategy::FULL;
    SolverType solver_type = SolverType::CG;
    PreconditionerType precond = PreconditionerType::JACOBI;
    
    // Parse arguments
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--mesh") == 0 && i + 1 < argc) {
            mesh_file = argv[++i];
        } else if (strcmp(argv[i], "--order") == 0 && i + 1 < argc) {
            order = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--refine") == 0 && i + 1 < argc) {
            refine_levels = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--assembly") == 0 && i + 1 < argc) {
            std::string asm_str = argv[++i];
            if (asm_str == "full") asm_strategy = AssemblyStrategy::FULL;
            else if (asm_str == "partial") asm_strategy = AssemblyStrategy::PARTIAL;
            else if (asm_str == "matfree") asm_strategy = AssemblyStrategy::MATRIX_FREE;
            else if (asm_str == "element") asm_strategy = AssemblyStrategy::ELEMENT_BY_ELEMENT;
        } else if (strcmp(argv[i], "--solver") == 0 && i + 1 < argc) {
            std::string solver_str = argv[++i];
            if (solver_str == "cg") {
                solver_type = SolverType::CG;
            } else if (solver_str == "gmres") {
                solver_type = SolverType::GMRES;
            } else if (solver_str == "amg") {
                solver_type = SolverType::CG;
                precond = PreconditionerType::AMG;  // AMG as preconditioner
            } else if (solver_str == "direct") {
                solver_type = SolverType::DIRECT;
            }
        } else if (strcmp(argv[i], "--precond") == 0 && i + 1 < argc) {
            std::string precond_str = argv[++i];
            if (precond_str == "none") precond = PreconditionerType::NONE;
            else if (precond_str == "jacobi") precond = PreconditionerType::JACOBI;
            else if (precond_str == "gs") precond = PreconditionerType::GAUSS_SEIDEL;
            else if (precond_str == "amg") precond = PreconditionerType::AMG;
        } else if (strcmp(argv[i], "--help") == 0) {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << "Options:\n"
                      << "  --mesh <file>       Mesh file (default: testdata/testmesh_cube.mesh)\n"
                      << "  --order <n>         FE order (default: 2)\n"
                      << "  --refine <n>        Refinement levels (default: 1)\n"
                      << "  --assembly <type>   Assembly strategy: full|partial|matfree|element (default: full)\n"
                      << "  --solver <type>     Solver type: cg|gmres|amg|direct (default: cg)\n"
                      << "  --precond <type>    Preconditioner: none|jacobi|gs|amg (default: jacobi)\n"
                      << "  --help              Show this help message\n";
            return 0;
        }
    }
    
    print_plugin_info();
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "Plugin-Based Poisson Solver" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Configuration:" << std::endl;
    std::cout << "  Mesh: " << mesh_file << std::endl;
    std::cout << "  Order: " << order << std::endl;
    std::cout << "  Refinement levels: " << refine_levels << std::endl;
    
    // Map enums to names for display
    auto get_asm_name = [](AssemblyStrategy s) {
        switch(s) {
            case AssemblyStrategy::FULL: return "Full Assembly";
            case AssemblyStrategy::PARTIAL: return "Partial Assembly";
            case AssemblyStrategy::MATRIX_FREE: return "Matrix-Free";
            case AssemblyStrategy::ELEMENT_BY_ELEMENT: return "Element-by-Element";
            default: return "Unknown";
        }
    };
    
    auto get_solver_name = [](SolverType s, PreconditionerType p) {
        std::string name;
        switch(s) {
            case SolverType::CG: name = "Conjugate Gradient"; break;
            case SolverType::GMRES: name = "GMRES"; break;
            case SolverType::DIRECT: name = "Direct (UMFPACK)"; break;
            default: name = "Unknown"; break;
        }
        if (p == PreconditionerType::AMG) {
            name += " + AMG";
        }
        return name;
    };
    
    auto get_precond_name = [](PreconditionerType p) {
        switch(p) {
            case PreconditionerType::NONE: return "None";
            case PreconditionerType::JACOBI: return "Jacobi";
            case PreconditionerType::GAUSS_SEIDEL: return "Gauss-Seidel";
            case PreconditionerType::AMG: return "AMG";
            case PreconditionerType::ILU: return "ILU";
            default: return "Unknown";
        }
    };
    
    std::cout << "  Assembly: " << get_asm_name(asm_strategy) << std::endl;
    std::cout << "  Solver: " << get_solver_name(solver_type, precond) << std::endl;
    std::cout << "  Preconditioner: " << get_precond_name(precond) << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Load and refine mesh
    Mesh mesh(mesh_file);
    int dim = mesh.Dimension();
    
    for (int i = 0; i < refine_levels; i++) {
        mesh.UniformRefinement();
    }
    
    std::cout << "Mesh statistics:" << std::endl;
    std::cout << "  Dimension: " << dim << std::endl;
    std::cout << "  Elements: " << mesh.GetNE() << std::endl;
    std::cout << "  Boundary elements: " << mesh.GetNBE() << std::endl;
    
    // Create finite element space
    H1_FECollection fec(order, dim);
    FiniteElementSpace fes(&mesh, &fec);
    
    std::cout << "  DOFs: " << fes.GetTrueVSize() << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Define problem: -Δu = f with homogeneous Dirichlet BC
    std::cout << "Assembling system..." << std::endl;
    
    // Create assembly plugin
    auto assembly = create_assembly_plugin(asm_strategy);
    
    // Bilinear form (Laplacian)
    BilinearForm a(&fes);
    ConstantCoefficient one(1.0);
    a.AddDomainIntegrator(new DiffusionIntegrator(one));
    
    // Assemble using plugin
    assembly->assemble(a);
    
    // Linear form (RHS)
    LinearForm b(&fes);
    ConstantCoefficient f_coeff(1.0);
    b.AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
    b.Assemble();
    
    // Essential boundary conditions (homogeneous Dirichlet)
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 1;  // Mark all boundaries as essential
    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    
    // Solution vector
    GridFunction x(&fes);
    x = 0.0;
    
    // Form linear system
    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);
    
    std::cout << "System statistics:" << std::endl;
    std::cout << "  Size: " << A.Height() << " x " << A.Width() << std::endl;
    std::cout << "  Non-zeros: " << A.NumNonZeroElems() << std::endl;
    std::cout << "  Assembly time: " << std::fixed << std::setprecision(3) 
              << assembly->metrics().setup_time_ms << " ms" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Solve system using plugin
    std::cout << "Solving system..." << std::endl;
    
    auto solver = create_solver_plugin(solver_type, precond);
    solver->set_tolerance(1e-9, 1e-12);
    solver->set_max_iterations(1000);
    solver->set_print_level(1);
    
    solver->setup(A);
    solver->solve(B, X);
    
    // Recover solution
    a.RecoverFEMSolution(X, b, x);
    
    // Print solver metrics
    const auto& metrics = solver->metrics();
    std::cout << "\nSolver statistics:" << std::endl;
    std::cout << "  Setup time: " << std::fixed << std::setprecision(3) 
              << metrics.setup_time_ms << " ms" << std::endl;
    std::cout << "  Solve time: " << std::fixed << std::setprecision(3) 
              << metrics.solve_time_ms << " ms" << std::endl;
    std::cout << "  Total time: " << std::fixed << std::setprecision(3) 
              << metrics.total_time_ms << " ms" << std::endl;
    std::cout << "  Iterations: " << metrics.iterations << std::endl;
    std::cout << "  Final residual: " << std::scientific << std::setprecision(2) 
              << metrics.residual << std::endl;
    
    // Compute solution statistics
    double min_val = x.Min();
    double max_val = x.Max();
    double l2_norm = x.Norml2();
    
    std::cout << "\nSolution statistics:" << std::endl;
    std::cout << "  Min value: " << std::fixed << std::setprecision(6) << min_val << std::endl;
    std::cout << "  Max value: " << std::fixed << std::setprecision(6) << max_val << std::endl;
    std::cout << "  L2 norm: " << std::scientific << std::setprecision(6) << l2_norm << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    // Save solution
    std::ofstream mesh_ofs("results/ex5_plugin_demo.mesh");
    mesh_ofs.precision(8);
    mesh.Print(mesh_ofs);
    
    std::ofstream sol_ofs("results/ex5_plugin_demo.gf");
    sol_ofs.precision(8);
    x.Save(sol_ofs);
    
    std::cout << "Solution saved to results/ex5_plugin_demo.{mesh,gf}" << std::endl;
    std::cout << "========================================\n" << std::endl;
    
    return 0;
}
