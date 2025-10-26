#include "mfem.hpp"
#include <fstream>
#include <iostream>

// Use the MFEM namespace
using namespace mfem;

int main(int argc, char *argv[])
{
    int myid = 0;

    // 2. Define the mesh and problem parameters
    const char *mesh_file = "testdata/testmesh1.msh";
    // Parse CmdLine arg (filename) is the user passes it
    if (argc > 1) {
        mesh_file = argv[1];
    }
    int order = 1; // P1 (linear) finite elements

    // 3. Read the mesh from the Gmsh file
    Mesh mesh(mesh_file);
    int dim = mesh.Dimension();
    std ::cout << "Mesh dimension: " << dim << "\n";
    if (myid == 0) {
        std::cout << "Mesh loaded: " << mesh.GetNE() << " elements.\n";
    }

    // 4. Define the Finite Element space (H1 for thermal)
    H1_FECollection fec(order, dim);
    FiniteElementSpace fes(&mesh, &fec);

    // 5. Define the Bilinear Form 'a(u,v)' (the left-hand side)
    // This is the integral of (k * grad(u) * grad(v))
    BilinearForm a(&fes);
    ConstantCoefficient k(401.0); // Thermal conductivity k=1
    a.AddDomainIntegrator(new DiffusionIntegrator(k));
    a.Assemble();

    // 6. Define the Linear Form 'b(v)' (the right-hand side)
    // This includes the heat source f
    LinearForm b(&fes);

    // Define a Gaussian heat source at the center
    Vector center(dim);
    Vector min_bb(dim), max_bb(dim);
    mesh.GetBoundingBox(min_bb, max_bb);
    for (int i = 0; i < dim; i++) {
        center(i) = (min_bb(i) + max_bb(i)) * 0.5;
    }
    
    // A C++ lambda function for the heat source f(x) returning a double
    auto gaussian_source = [&](const Vector &x) -> double {
        Vector x_minus_c(x);
        x_minus_c -= center;
        double r = x_minus_c.Norml2();
        double r2 = r * r;
        double A = 1e8; // Source amplitude
        double b = 50.0; // Source "width"
        return A * exp(-b * r2);
    };
    FunctionCoefficient f_coeff(gaussian_source);
    
    b.AddDomainIntegrator(new DomainLFIntegrator(f_coeff));
    b.Assemble();

    // 7. Define the solution GridFunction 'x' (Temperature)
    GridFunction x(&fes);
    x = 0.0; // Initialize temperature to 0

    // 8. Set up the Dirichlet Boundary Conditions
    // We assume:
    // - Boundary Attribute 1 = "hot" side (T=1.0)
    // - Boundary Attribute 2 = "cold" side (T=0.0)

    // Get a list of all essential (Dirichlet) degrees of freedom
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[1-1] = 1; // Attribute 1
    ess_bdr[2-1] = 1; // Attribute 2

    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // Apply the T=1.0 and T=0.0 values
    ConstantCoefficient hot_temp(100.0);
    ConstantCoefficient cold_temp(0.0);
    
    for (auto &i : mesh.bdr_attributes) {
        std::cout << "Boundary attribute: " << i << "\n";
    }

    Array<int> hot_bdr(mesh.bdr_attributes.Max());
    hot_bdr = 0; hot_bdr[1-1] = 1; // Set only Attr 1
    x.ProjectBdrCoefficient(hot_temp, hot_bdr);

    Array<int> cold_bdr(mesh.bdr_attributes.Max());
    cold_bdr = 0; cold_bdr[2-1] = 1; // Set only Attr 2
    x.ProjectBdrCoefficient(cold_temp, cold_bdr);

    // 9. Form the final linear system A*X = B
    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, b, A, X, B);

    // 10. Solve the linear system A*X = B
    // We'll use a simple preconditioned conjugate gradient (PCG) solver
    if (myid == 0) {
        std::cout << "Solving linear system...\n";
    }
    GSSmoother M(A); // Use a simple smoother as a preconditioner
    // Use the MFEM CG solver
    CGSolver solver;
    solver.SetOperator(A);
    solver.SetRelTol(1e-12);
    solver.SetMaxIter(2000);
    solver.SetPrintLevel(3); // Print solver progress
    solver.SetPreconditioner(M);
    solver.Mult(B, X);

    // 11. Recover the solution into the GridFunction 'x'
    a.RecoverFEMSolution(X, b, x);

    // 12. Save the solution for ParaView
    if (myid == 0) {
        std::cout << "Saving solution to ParaView files...\n";
    }
    ParaViewDataCollection paraview_dc("MyThermalSolution", &mesh);
    paraview_dc.SetLevelsOfDetail(order); // Save high-order data
    paraview_dc.RegisterField("Temperature", &x);
    paraview_dc.SetDataFormat(mfem::VTKFormat::ASCII);
    paraview_dc.Save();

    // 13. Finalize
    if (myid == 0) {
        std::cout << "Done. Open 'MyThermalSolution/MyThermalSolution.pvd' in ParaView.\n";
    }
    return 0;
}