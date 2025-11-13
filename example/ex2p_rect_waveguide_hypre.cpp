// Adapted from ex2_rect_waveguide.cpp
// Uses HypreAmgSolver (hpcfem wrapper) instead of GMRES

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <complex>
#include <vector>
#include <cmath>

// Hypre AMG solver wrapper
#include "hpcfem/solvers/solver_hypre_amg.hpp"

// Rectangular waveguide example (frequency-domain, single-frequency excitation)
// - Uses Nedelec (ND) finite elements (H(curl))
// - Ports: boundary attribute 3 (input) and 4 (output/termination)
// - PEC on attributes 1,2,5,6 (tangential E = 0)
// - Modal excitation: TE10 on port 1 (simple analytical profile)
// - Mesh: testdata/testmesh_bar.mesh
// - Output: results/ex8_rect_waveguide_hypre

using namespace mfem;

template <typename T> T pow2(const T &x) { return x*x; }

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

// (PML/helpers copied from ex2; unchanged)
class PML
{
private:
    Mesh *mesh;
    int dim;
    Array2D<real_t> length; // PML length per direction (left/right)
    Array2D<real_t> comp_dom_bdr;
    Array2D<real_t> dom_bdr;
    Array<int> elems; // 0: in PML, 1: computational domain

    void SetBoundaries();

public:
    PML(Mesh *mesh_, Array2D<real_t> length_);
    Array2D<real_t> GetCompDomainBdr() { return comp_dom_bdr; }
    Array2D<real_t> GetDomainBdr() { return dom_bdr; }
    Array<int> * GetMarkedPMLElements() { return &elems; }

    void SetAttributes(Mesh *mesh_);
    void StretchFunction(const Vector &x, std::vector<std::complex<real_t>> &dxs);
};

class PMLDiagMatrixCoefficient : public VectorCoefficient
{
private:
    PML * pml = nullptr;
    void (*Function)(const Vector &, PML *, Vector &);
public:
    PMLDiagMatrixCoefficient(int dim, void(*F)(const Vector &, PML *, Vector &),
                             PML * pml_)
        : VectorCoefficient(dim), pml(pml_), Function(F)
    {}

    using VectorCoefficient::Eval;

    void Eval(Vector &K, ElementTransformation &T, const IntegrationPoint &ip) override
    {
        real_t x[3];
        Vector transip(x, 3);
        T.Transform(ip, transip);
        K.SetSize(vdim);
        (*Function)(transip, pml, K);
    }
};

void detJ_JT_J_inv_abs(const Vector &x, PML * pml, Vector &D);
void detJ_inv_JT_J_abs(const Vector &x, PML * pml, Vector &D);

PML::PML(Mesh *mesh_, Array2D<real_t> length_)
    : mesh(mesh_), length(length_)
{
    dim = mesh->Dimension();
    SetBoundaries();
}

void PML::SetBoundaries()
{
    comp_dom_bdr.SetSize(dim, 2);
    dom_bdr.SetSize(dim, 2);
    Vector pmin, pmax;
    mesh->GetBoundingBox(pmin, pmax);
    for (int i = 0; i < dim; i++)
    {
        dom_bdr(i, 0) = pmin(i);
        dom_bdr(i, 1) = pmax(i);
        comp_dom_bdr(i, 0) = dom_bdr(i, 0) + length(i, 0);
        comp_dom_bdr(i, 1) = dom_bdr(i, 1) - length(i, 1);
    }
}

void PML::SetAttributes(Mesh *mesh_)
{
    int nrelem = mesh_->GetNE();
    elems.SetSize(nrelem);

    for (int i = 0; i < nrelem; ++i)
    {
        elems[i] = 1; // default: computational domain
        bool in_pml = false;
        Element *el = mesh_->GetElement(i);
        Array<int> vertices;
        el->GetVertices(vertices);
        int nrvert = vertices.Size();

        // Initialize element attribute to 1
        el->SetAttribute(1);

        for (int iv = 0; iv < nrvert; ++iv)
        {
            int vert_idx = vertices[iv];
            real_t *coords = mesh_->GetVertex(vert_idx);
            for (int comp = 0; comp < dim; ++comp)
            {
                if (coords[comp] > comp_dom_bdr(comp, 1) ||
                    coords[comp] < comp_dom_bdr(comp, 0))
                {
                    in_pml = true;
                    break;
                }
            }
            if (in_pml) { break; }
        }
        if (in_pml)
        {
            elems[i] = 0;
            el->SetAttribute(2);
        }
    }
    mesh_->SetAttributes();
}

void PML::StretchFunction(const Vector &x, std::vector<std::complex<real_t>> &dxs)
{
    constexpr std::complex<real_t> zi = std::complex<real_t>(0., 1.);
    real_t n = 2.0;
    real_t c = 5.0;
    real_t coeff;

    for (int i = 0; i < dim; ++i)
    {
        dxs[i] = std::complex<real_t>(1.0, 0.0);
        if (x(i) >= comp_dom_bdr(i, 1))
        {
            coeff = n * c / pow(length(i, 1), n);
            dxs[i] = std::complex<real_t>(1_r, 0.0) + zi * coeff * std::abs(pow(x(i) - comp_dom_bdr(i, 1), n - 1_r));
        }
        if (x(i) <= comp_dom_bdr(i, 0))
        {
            coeff = n * c / pow(length(i, 0), n);
            dxs[i] = std::complex<real_t>(1_r, 0.0) + zi * coeff * std::abs(pow(x(i) - comp_dom_bdr(i, 0), n - 1_r));
        }
    }
}

void detJ_JT_J_inv_abs(const Vector &x, PML * pml, Vector &D)
{
    int pdim = D.Size();
    std::vector<std::complex<real_t>> dxs(pdim);
    std::complex<real_t> det = std::complex<real_t>(1.0, 0.0);
    pml->StretchFunction(x, dxs);
    for (int i = 0; i < pdim; ++i) { det *= dxs[i]; }
    for (int i = 0; i < pdim; ++i)
    {
        D(i) = std::abs(det / pow2(dxs[i]));
    }
}

void detJ_inv_JT_J_abs(const Vector &x, PML * pml, Vector &D)
{
    int pdim = D.Size();
    std::vector<std::complex<real_t>> dxs(pdim);
    std::complex<real_t> det = std::complex<real_t>(1.0, 0.0);
    pml->StretchFunction(x, dxs);
    for (int i = 0; i < pdim; ++i) { det *= dxs[i]; }
    if (pdim == 2)
    {
        D = std::abs(std::complex<real_t>(1_r, 0.0) / det);
    }
    else
    {
        for (int i = 0; i < pdim; ++i)
        {
            D(i) = std::abs(pow2(dxs[i]) / det);
        }
    }
}

int main(int argc, char *argv[])
{
    const char *mesh_file = "testdata/testmesh_bar.mesh";
    const char *output_dir = "results/ex2p_rect_waveguide_hypre";
    if (argc > 1)
    {
        mesh_file = argv[1];
    }
    if (argc > 2)
    {
        output_dir = argv[2];
    }

    // Physical & frequency parameters
    const double a = 0.02286;      // x length (m)
    const double b = 0.01016;      // y length (m)
    const double freq = 10e9;      // 10 GHz
    const double c0 = 299792458.0; // speed of light
    const double omega = 2.0 * M_PI * freq;
    const double k0 = omega / c0; // free-space wavenumber
    const double eps_r = 1.0;     // relative permittivity (vacuum)
    const double mu_r = 1.0;      // relative permeability

    // Finite element order
    int order = 1;

#ifdef MFEM_USE_MPI
    MPI_Session mpi(argc, argv);
    int myid = mpi.WorldRank();
#else
    int myid = 0;
#endif

    // Read mesh (serial first, then convert to ParMesh if MPI enabled)
    Mesh *mesh = new Mesh(mesh_file);
    int dim = mesh->Dimension();
    if (dim != 3)
    {
        std::cout << "ex8_rect_waveguide_hypre requires a 3D mesh. Mesh dimension: " << dim << std::endl;
        return 1;
    }

#ifdef MFEM_USE_MPI
    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
#else
    Mesh *pmesh = mesh;
#endif

    if (myid == 0)
    {
        std::cout << "Mesh loaded: " << pmesh->GetNE() << " elements." << std::endl;
    }

    // -------------------- Setup PML region --------------------
    Array2D<real_t> length(dim, 2);
    length = 0.0;
    if (dim >= 1)
    {
        length(0, 1) = 0.001; // 2 mm PML thickness (tunable)
    }
    // PML accepts Mesh* and ParMesh derives from Mesh, so pass pmesh
    PML pml(pmesh, length);
    pml.SetAttributes(pmesh);

#ifdef MFEM_USE_MPI
    ND_FECollection fec(order, dim);
    ParFiniteElementSpace fespace(pmesh, &fec);
    HYPRE_BigInt true_size = fespace.GlobalTrueVSize();
    if (myid == 0) std::cout << "Number of H(curl) dofs: " << true_size << std::endl;

    ParBilinearForm aform(&fespace);
    ParLinearForm rhs(&fespace);
    ParGridFunction x(&fespace);
#else
    ND_FECollection fec(order, dim);
    FiniteElementSpace fespace(pmesh, &fec);
    int true_size = fespace.GetTrueVSize();
    std::cout << "Number of H(curl) dofs: " << true_size << std::endl;

    BilinearForm aform(&fespace);
    LinearForm rhs(&fespace);
    GridFunction x(&fespace);
#endif

    // Bilinear form: curl-curl - k0^2 * eps * I
    Array<int> attr(pmesh->attributes.Max());
    Array<int> attrPML(pmesh->attributes.Max());
    attr = 0; attr[0] = 1;
    attrPML = 0;
    if (pmesh->attributes.Max() > 1)
    {
        attrPML[1] = 1; // attribute 2 -> PML region
    }

    ConstantCoefficient muinv(1.0 / mu_r);
    double mass_coef_val = -(k0 * k0 * eps_r);
    ConstantCoefficient omeg(mass_coef_val);

    RestrictedCoefficient restr_muinv(muinv, attr);
    RestrictedCoefficient restr_omeg(omeg, attr);

    aform.AddDomainIntegrator(new CurlCurlIntegrator(restr_muinv));
    aform.AddDomainIntegrator(new VectorFEMassIntegrator(restr_omeg));

    int cdim = (dim == 2) ? 1 : dim;
    PMLDiagMatrixCoefficient pml_c1_abs(cdim, detJ_inv_JT_J_abs, &pml);
    ScalarVectorProductCoefficient c1_abs(muinv, pml_c1_abs);
    VectorRestrictedCoefficient restr_c1_abs(c1_abs, attrPML);

    ConstantCoefficient absomeg(std::abs(mass_coef_val));
    PMLDiagMatrixCoefficient pml_c2_abs(dim, detJ_JT_J_inv_abs, &pml);
    ScalarVectorProductCoefficient c2_abs(absomeg, pml_c2_abs);
    VectorRestrictedCoefficient restr_c2_abs(c2_abs, attrPML);

    aform.AddDomainIntegrator(new CurlCurlIntegrator(restr_c1_abs));
    aform.AddDomainIntegrator(new VectorFEMassIntegrator(restr_c2_abs));

    // Solution container
    x = 0.0;

    Array<int> ess_bdr(pmesh->bdr_attributes.Max());
    ess_bdr = 0;
    ess_bdr[1 - 1] = 1;
    ess_bdr[2 - 1] = 1;
    ess_bdr[5 - 1] = 1;
    ess_bdr[6 - 1] = 1;
    ess_bdr[3 - 1] = 1; // port 3: driven

    Array<int> ess_tdof_list;
#ifdef MFEM_USE_MPI
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
#else
    fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
#endif

    double amplitude = 1.0;
    TE10Modal te10(a, b, amplitude);
    Array<int> port3(pmesh->bdr_attributes.Max());
    port3 = 0;
    port3[3 - 1] = 1;
    x.ProjectBdrCoefficientTangent(te10, port3);

    // Assemble
    aform.Assemble();

    std::cout << "Projected boundary x norm (L2): " << x.Norml2()
              << ", Linf: " << x.Normlinf() << std::endl;

    // Form linear system and solve
#ifdef MFEM_USE_MPI
    HypreParMatrix A;
    Vector B, X;
    rhs = 0.0;
    rhs.Assemble();
    aform.FormLinearSystem(ess_tdof_list, x, rhs, A, X, B);
    if (myid == 0)
    {
        std::cout << "Assembled parallel matrix. True DOFs: " << A.GetGlobalNumRows() << std::endl;
    }

    std::cout << "Solving linear system with HypreAmgSolver..." << std::endl;
    hpcfem::HypreAmgSolver amg_solver(1e-8, 2000, 0);
    amg_solver.solve(A, B, X);
    if (myid == 0)
    {
        std::cout << "AMG iterations: " << amg_solver.getNumIterations() << std::endl;
        std::cout << "AMG final norm: " << amg_solver.getFinalNorm() << std::endl;
    }

    aform.RecoverFEMSolution(X, rhs, x);

    // Save to ParaView (parallel)
    ParaViewDataCollection paraview_dc(output_dir, pmesh);
    paraview_dc.SetLevelsOfDetail(1);
    paraview_dc.RegisterField("E_field", &x);
    paraview_dc.SetDataFormat(mfem::VTKFormat::ASCII);
    paraview_dc.Save();
#else
    SparseMatrix A;
    Vector B, X;
    rhs = 0.0;
    rhs.Assemble();
    aform.FormLinearSystem(ess_tdof_list, x, rhs, A, X, B);

    std::cout << "Assembled matrix size: " << A.Height() << " x " << A.Width() << std::endl;
    std::cout << "RHS (reduced) norm L2: " << B.Norml2() << ", Linf: " << B.Normlinf() << std::endl;

    std::cout << "Solving linear system with HypreAmgSolver..." << std::endl;
    hpcfem::HypreAmgSolver amg_solver(1e-8, 2000, 0);
    amg_solver.solve(A, B, X);
    std::cout << "AMG iterations: " << amg_solver.getNumIterations() << std::endl;
    std::cout << "AMG final norm: " << amg_solver.getFinalNorm() << std::endl;

    // Recover solution
    aform.RecoverFEMSolution(X, rhs, x);

    // Save to ParaView
    std::cout << "Saving solution to ParaView files..." << std::endl;
    ParaViewDataCollection paraview_dc(output_dir, pmesh);
    paraview_dc.SetLevelsOfDetail(1);
    paraview_dc.RegisterField("E_field", &x);
    paraview_dc.SetDataFormat(mfem::VTKFormat::ASCII);
    paraview_dc.Save();
#endif

    if (myid == 0)
    {
        std::cout << "Done. Output in '" << output_dir << "'." << std::endl;
    }

    return 0;
}
