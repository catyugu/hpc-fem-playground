#include "mfem.hpp"
#include <fstream>
#include <iostream>

// Example: steady-state scalar diffusion/heat equation with mixed BCs
// Problem summary (what, where):
// - Domain: mesh (Gmsh-format) provided by `mesh_file`.
// - Unknown: scalar temperature T in H^1 (continuous finite elements).
// - PDE (strong form): -div(k grad T) = 0 in the domain.
// - Boundary conditions:
//     * Attribute 1 : Dirichlet T = T_cold (essential)
//     * Attribute 6 : Robin (convective) -k dT/dn = h*(T - T_inf)
//     * Attribute 2 : Neumann prescribed inward flux q (natural)
// Implementation outline (high level):
//  1) Read mesh and build H1 FE space (order = 1 here).
//  2) Build BilinearForm `a` with DiffusionIntegrator(k).
//     - Add MassIntegrator(h) on the Robin boundary to include h*u*v term.
//  3) Build LinearForm `lf` for RHS contributions:
//     - Robin contribution: 
//           \int_{Gamma_robin} h * T_inf * v
//       implemented with `BoundaryLFIntegrator(h*T_inf)` on the Robin markers.
//     - Neumann contribution: 
//           \int_{Gamma_neu} q * v
//       implemented with `BoundaryLFIntegrator(q)` on the Neumann markers.
//  4) Apply Dirichlet BCs by setting essential boundary attributes and
//     projecting their values into `x` (initial GridFunction / known part).
//  5) Call `a.Assemble()` and `lf.Assemble()` then `a.FormLinearSystem(ess_tdof_list, x, lf, A, X, B)`.
//  6) Solve the linear system (CG + preconditioner) and recover FEM solution
//     with `a.RecoverFEMSolution(X, lf, x)`.
//  7) Save the `x` GridFunction via `ParaViewDataCollection`.

// Use the MFEM namespace
using namespace mfem;

int main(int argc, char *argv[])
{
    int myid = 0;

    // --- (I) Inputs and user options -------------------------------------------------
    const char *mesh_file = "testdata/testmesh_cube.mesh";
    const char *output_dir = "results/ex1_cube_mixed_bc";
    if (argc > 1) { mesh_file = argv[1]; }
    if (argc > 2) { output_dir = argv[2]; }
    int order = 1; // finite element order (P1 = linear)

    // --- (II) Mesh and FE space -----------------------------------------------------
    Mesh mesh(mesh_file);
    int dim = mesh.Dimension();
    std::cout << "Mesh dimension: " << dim << "\n";
    if (myid == 0) { std::cout << "Mesh loaded: " << mesh.GetNE() << " elements." << "\n"; }

    H1_FECollection fec(order, dim);
    FiniteElementSpace fes(&mesh, &fec);

    // --- (III) Bilinear form (left-hand side): diffusion + Robin boundary mass term
    BilinearForm a(&fes);
    ConstantCoefficient k(401.0); // thermal conductivity
    a.AddDomainIntegrator(new DiffusionIntegrator(k));

    // --- (IV) Solution container / initial guess -----------------------------------
    GridFunction x(&fes);
    x = 0.0; // initialize (will contain projected Dirichlet values)

    // --- (V) Essential (Dirichlet) BCs ---------------------------------------------
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[1-1] = 1; // Attribute 1 is essential (cold Dirichlet)
    // Note: attribute 6 is NOT made essential here: it is handled as Robin (natural part)
    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    ConstantCoefficient cold_temp(293.15);
    ConstantCoefficient hot_temp(393.15); // ambient/hot world temperature used for projection/info

    // Print available boundary attributes for clarity
    for (auto &i : mesh.bdr_attributes) { std::cout << "Boundary attribute: " << i << "\n"; }

    Array<int> cold_bdr(mesh.bdr_attributes.Max());
    cold_bdr = 0; cold_bdr[1-1] = 1; // attribute 1 only
    x.ProjectBdrCoefficient(cold_temp, cold_bdr);

    Array<int> hot_bdr(mesh.bdr_attributes.Max());
    hot_bdr = 0; hot_bdr[6-1] = 1; // attribute 6 (Robin)

    // Robin BC setup: add h * u * v to the bilinear form on the robin marker
    double h_coeff_val = 100.0; // heat transfer coefficient (W/m^2/K)
    ConstantCoefficient h_coeff(h_coeff_val);
    a.AddBoundaryIntegrator(new MassIntegrator(h_coeff), hot_bdr);

    // --- (VI) Linear form (RHS): Robin and Neumann contributions --------------------
    // Robin RHS: h * T_inf
    ConstantCoefficient robin_rhs_coef(h_coeff_val * 393.15);
    LinearForm lf(&fes);
    lf.AddBoundaryIntegrator(new BoundaryLFIntegrator(robin_rhs_coef), hot_bdr);

    // Neumann BC on attribute 2: prescribed inward heat flux q (W/m^2)
    // Adds: \int_{Gamma_neu} q * v  (positive q -> heat entering the domain)
    Array<int> neumann_bdr(mesh.bdr_attributes.Max());
    neumann_bdr = 0; neumann_bdr[2-1] = 1; // attribute 2
    double inward_flux_val = 5e4; // W/m^2
    ConstantCoefficient neumann_flux(inward_flux_val);
    lf.AddBoundaryIntegrator(new BoundaryLFIntegrator(neumann_flux), neumann_bdr);

    // Assemble system contributions
    a.Assemble();
    lf.Assemble();

    // --- (VII) Form linear system and solve ----------------------------------------
    SparseMatrix A; Vector B, X;
    a.FormLinearSystem(ess_tdof_list, x, lf, A, X, B);

    if (myid == 0) { std::cout << "Solving linear system...\n"; }
    GSSmoother M(A);
    CGSolver solver;
    solver.SetOperator(A);
    solver.SetRelTol(1e-12);
    solver.SetMaxIter(2000);
    solver.SetPrintLevel(3);
    solver.SetPreconditioner(M);
    solver.Mult(B, X);

    // Recover the full FEM solution into x
    a.RecoverFEMSolution(X, lf, x);

    // --- (VIII) Save ---------------------------------------------------------------
    if (myid == 0) { std::cout << "Saving solution to ParaView files...\n"; }
    ParaViewDataCollection paraview_dc(output_dir, &mesh);
    paraview_dc.SetLevelsOfDetail(1);
    paraview_dc.RegisterField("Temperature", &x);
    paraview_dc.SetDataFormat(mfem::VTKFormat::ASCII);
    paraview_dc.Save();

    if (myid == 0) { std::cout << "Done. Open '" << output_dir << "/" << output_dir << ".pvd' in ParaView.\n"; }
    return 0;
}