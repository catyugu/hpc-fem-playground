#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

// Rectangular waveguide example (frequency-domain, single-frequency excitation)
// - Uses Nedelec (ND) finite elements (H(curl))
// - Ports: boundary attribute 3 (input) and 4 (output/termination)
// - PEC on attributes 1,2,5,6 (tangential E = 0)
// - Modal excitation: TE10 on port 1 (simple analytical profile)
// - Mesh: testdata/testmesh_bar.mesh
// - Output: results/ex2_rect_waveguide

using namespace mfem;

class TE10Modal : public VectorCoefficient
{
private:
    double a_; // width (x)
    double b_; // height (y)
    double amp_;

public:
    // VectorCoefficient expects the vector dimension (3 for 3D)
    TE10Modal(double a, double b, double amp)
        : VectorCoefficient(3), a_(a), b_(b), amp_(amp) {}

    // Evaluate using physical coordinates via the ElementTransformation
    void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip)
    {
        Vector X(3);
        T.Transform(ip, X);
        double x = X(0);
        double y = X(1);
        V.SetSize(3);
        // TE10 dominant transverse-E field (approx.): E = (0, E_y, 0)
        // E_y ~ sin(pi * x / a)
        double Ey = std::sin(M_PI * x / a_);
        V(0) = 0.0;       // E_x
        V(1) = amp_ * Ey; // E_y
        V(2) = 0.0;       // E_z
    }
};

int main(int argc, char *argv[])
{
    int myid = 0;

    const char *mesh_file = "testdata/testmesh_bar.mesh";
    const char *output_dir = "results/ex2_rect_waveguide";
    if (argc > 1)
    {
        mesh_file = argv[1];
    }
    if (argc > 2)
    {
        output_dir = argv[2];
    }

    // Physical & frequency parameters
    const double a = 0.01016;      // x length (m)
    const double b = 0.02286;      // y length (m)
    const double freq = 10e9;      // 10 GHz
    const double c0 = 299792458.0; // speed of light
    const double omega = 2.0 * M_PI * freq;
    const double k0 = omega / c0; // free-space wavenumber
    const double eps_r = 1.0;     // relative permittivity (vacuum)
    const double mu_r = 1.0;      // relative permeability

    // Finite element order
    int order = 1;

    // Read mesh
    Mesh mesh(mesh_file);
    int dim = mesh.Dimension();
    if (dim != 3)
    {
        std::cout << "ex2_rect_waveguide requires a 3D mesh. Mesh dimension: " << dim << std::endl;
        return 1;
    }

    if (myid == 0)
    {
        std::cout << "Mesh loaded: " << mesh.GetNE() << " elements." << std::endl;
    }

    // Debug: report boundary attribute counts
    int nbdr = mesh.GetNBE();
    Array<int> bdr_attr_counts(mesh.bdr_attributes.Max());
    bdr_attr_counts = 0;
    for (int be = 0; be < nbdr; be++)
    {
        Element *b = mesh.GetBdrElement(be);
        int a = b->GetAttribute();
        if (a >= 1 && a <= bdr_attr_counts.Size())
        {
            bdr_attr_counts[a - 1]++;
        }
    }
    std::cout << "Boundary elements: " << nbdr << std::endl;
    for (int i = 0; i < bdr_attr_counts.Size(); i++)
    {
        if (bdr_attr_counts[i] > 0)
        {
            std::cout << "  attr " << (i + 1) << ": " << bdr_attr_counts[i] << " faces" << std::endl;
        }
    }

    // ND (Nedelec) finite element space for H(curl)
    ND_FECollection fec(order, dim);
    FiniteElementSpace fes(&mesh, &fec);
    if (myid == 0)
    {
        std::cout << "Number of H(curl) dofs: " << fes.GetTrueVSize() << std::endl;
    }

    // Bilinear form: curl-curl - k0^2 * eps * I
    BilinearForm aform(&fes);
    ConstantCoefficient one_over_mu(1.0 / mu_r);
    aform.AddDomainIntegrator(new CurlCurlIntegrator(one_over_mu));
    // mass term with negative sign to form (curl curl - k^2 eps) weak form
    double mass_coef_val = -(k0 * k0 * eps_r);
    ConstantCoefficient mass_coef(mass_coef_val);
    aform.AddDomainIntegrator(new VectorFEMassIntegrator(mass_coef));

    // Solution container
    GridFunction x(&fes);
    x = 0.0;

    // Essential (Dirichlet) BCs:
    // - PEC on attributes 1,2,5,6 (tangential E = 0)
    // - Ports on attributes 3 (drive) and 4 (termination)
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 0;
    // PEC (set to zero tangential E)
    ess_bdr[1 - 1] = 1;
    ess_bdr[2 - 1] = 1;
    ess_bdr[5 - 1] = 1;
    ess_bdr[6 - 1] = 1;
    // Ports (will be projected below)
    ess_bdr[3 - 1] = 1; // port 3: driven
    ess_bdr[4 - 1] = 1; // port 4: termination

    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // Project port fields
    // Drive amplitude (arbitrary scale)
    double amplitude = 1.0;
    TE10Modal te10(a, b, amplitude);
    Array<int> port1(mesh.bdr_attributes.Max());
    port1 = 0;
    port1[1 - 1] = 1;
    // For H(curl) (ND) spaces, project the tangential component on the boundary
    // Drive port: attribute 3
    Array<int> port3(mesh.bdr_attributes.Max());
    port3 = 0;
    port3[3 - 1] = 1;
    x.ProjectBdrCoefficientTangent(te10, port3);

    // Port 6 set to zero (termination)
    // Termination port: attribute 4 (set to zero tangential field)
    Array<int> port4(mesh.bdr_attributes.Max());
    port4 = 0;
    port4[4 - 1] = 1;
    // Create a zero vector coefficient for termination
    Vector zvec(3);
    zvec = 0.0;
    VectorConstantCoefficient zero_vec(zvec);
    x.ProjectBdrCoefficientTangent(zero_vec, port4);

    // Assemble
    aform.Assemble();

    // Debug: report norm of the projected boundary GridFunction (should be non-zero)
    std::cout << "Projected boundary x norm (L2): " << x.Norml2()
              << ", Linf: " << x.Normlinf() << std::endl;

    // Form linear system and solve (indefinite Helmholtz-like). We'll use GMRES.
    SparseMatrix A;
    Vector B, X;
    LinearForm rhs(&fes);
    rhs = 0.0;
    rhs.Assemble();
    aform.FormLinearSystem(ess_tdof_list, x, rhs, A, X, B);

    std::cout << "Assembled matrix size: " << A.Height() << " x " << A.Width() << std::endl;
    std::cout << "RHS (reduced) norm L2: " << B.Norml2() << ", Linf: " << B.Normlinf() << std::endl;

    // Try a direct solver if available at compile time; otherwise fall back to GMRES
    bool solved = false;
#if defined(MFEM_USE_SUITESPARSE)
    if (myid == 0)
    {
        std::cout << "Solving linear system with UMFPACK (direct) ..." << std::endl;
    }
    UMFPackSolver umf;
    umf.SetOperator(A);
    umf.Mult(B, X);
    solved = true;
#elif defined(MFEM_USE_MUMPS)
    if (myid == 0)
    {
        std::cout << "Solving linear system with MUMPS (direct) ..." << std::endl;
    }
    // MUMPSSolver requires MPI; construct with MPI_COMM_WORLD
#ifdef MFEM_USE_MPI
    MUMPSSolver mumps(MPI_COMM_WORLD);
    mumps.SetOperator(A);
    mumps.Mult(B, X);
    solved = true;
#endif
#elif defined(MFEM_USE_SUPERLU)
    if (myid == 0)
    {
        std::cout << "Solving linear system with SuperLU_DIST (direct) ..." << std::endl;
    }
#ifdef MFEM_USE_MPI
    SuperLUSolver slu(MPI_COMM_WORLD);
    slu.SetOperator(A);
    slu.Mult(B, X);
    solved = true;
#endif
#endif

    if (!solved)
    {
        if (myid == 0)
        {
            std::cout << "Solving linear system (GMRES)..." << std::endl;
        }
        GSSmoother M(A);
        GMRESSolver solver;
        solver.SetOperator(A);
        solver.SetRelTol(1e-8);
        solver.SetAbsTol(1e-6);
        solver.SetMaxIter(500);
        solver.SetPrintLevel(1);
        solver.SetPreconditioner(M);
        solver.Mult(B, X);
    }

    // Recover solution
    aform.RecoverFEMSolution(X, rhs, x);

    // Save to ParaView
    if (myid == 0)
    {
        std::cout << "Saving solution to ParaView files..." << std::endl;
    }
    ParaViewDataCollection paraview_dc(output_dir, &mesh);
    paraview_dc.SetLevelsOfDetail(1);
    paraview_dc.RegisterField("E_field", &x);
    paraview_dc.SetDataFormat(mfem::VTKFormat::ASCII);
    paraview_dc.Save();

    if (myid == 0)
    {
        std::cout << "Done. Output in '" << output_dir << "'." << std::endl;
    }
    return 0;
}
