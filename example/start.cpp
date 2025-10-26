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

    // 7. Define the solution GridFunction 'x' (Temperature)
    GridFunction x(&fes);
    x = 0.0; // Initialize temperature to 0

    // 8. Set up the Dirichlet Boundary Conditions
    Array<int> ess_bdr(mesh.bdr_attributes.Max());

    ess_bdr = 0;
    ess_bdr[3-1] = 1; // Attribute 3
    ess_bdr[13-1] = 1; // Attribute 13

    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    ConstantCoefficient cold_temp(293.15);
    ConstantCoefficient hot_temp(1220.0);
    
    
    for (auto &i : mesh.bdr_attributes) {
        std::cout << "Boundary attribute: " << i << "\n";
    }

    Array<int> cold_bdr(mesh.bdr_attributes.Max());
    cold_bdr = 0; cold_bdr[3-1] = 1; // Set only Attr 3
    x.ProjectBdrCoefficient(cold_temp, cold_bdr);

    Array<int> hot_bdr(mesh.bdr_attributes.Max());
    hot_bdr = 0; hot_bdr[13-1] = 1; // Set only Attr 13
    x.ProjectBdrCoefficient(hot_temp, hot_bdr);

    // 9. Form the final linear system A*X = B
    SparseMatrix A;
    Vector B, X;
    // No load, so b is zero
    Vector b(fes.GetTrueVSize());
    b = 0.0;
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