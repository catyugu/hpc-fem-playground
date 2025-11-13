// Simple, clean 3D driven-frequency rectangular waveguide example.
// Solves: curl (mu^{-1} curl E) - k0^2 eps E = RHS (frequency domain)
// - Nedelec H(curl) elements, complex sesquilinear form
// - Port 3: driven TE10 (Dirichlet tangential), Port 4: modal termination

#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <complex>
#include <string>

using namespace std;
using namespace mfem;

// Simple TE10 modal profile (real-valued) on a rectangular port:
class TE10Modal : public VectorCoefficient
{
public:
    TE10Modal(double a, double amp = 1.0) : VectorCoefficient(3), a_(a), amp_(amp) {}
    void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip) override
    {
        Vector X(3); T.Transform(ip, X);
        const double x = X(0);
        V.SetSize(3);
        V(0) = 0.0;
        V(1) = amp_ * sin(M_PI * x / a_); // Ey shape of TE10
        V(2) = 0.0;
    }
private:
    double a_, amp_;
};

// Modal admittance matrix = Y * (emode * emode^T)
class ModalMatrixCoefficient : public MatrixCoefficient
{
public:
    ModalMatrixCoefficient(int dim, VectorCoefficient &mode, double Y)
        : MatrixCoefficient(dim), mode_(&mode), Y_(Y) {}
    void Eval(DenseMatrix &M, ElementTransformation &T, const IntegrationPoint &ip) override
    {
        Vector e(3); mode_->Eval(e, T, ip);
        M.SetSize(3,3);
        for (int i=0;i<3;i++) for (int j=0;j<3;j++) M(i,j) = Y_ * e(i) * e(j);
    }
private:
    VectorCoefficient *mode_;
    double Y_;
};

int main(int argc, char *argv[])
{
    const char *mesh_file = "testdata/testmesh_bar.mesh";
    const char *output_dir = "results/ex2_rect_waveguide";
    if (argc > 1) mesh_file = argv[1];
    if (argc > 2) output_dir = argv[2];

    // Physical parameters (10 GHz rectangular waveguide, WR-90-ish)
    const double a = 0.02286;      // width (x)
    const double b = 0.01016;      // height (y)
    const double freq = 10e9;
    const double c0 = 299792458.0;
    const double omega = 2.0 * M_PI * freq;
    const double k0 = omega / c0;
    const double eps_r = 1.0;
    const double mu_r = 1.0;

    // FE order (small for testing)
    int order = 1;

    // Read mesh
    Mesh mesh(mesh_file);
    if (mesh.Dimension() != 3)
    {
        cerr << "Error: mesh must be 3D." << endl;
        return 1;
    }

    // H(curl) space
    ND_FECollection fec(order, 3);
    FiniteElementSpace fes(&mesh, &fec);

    cout << "Mesh elements: " << mesh.GetNE()
         << ", H(curl) true DOFs: " << fes.GetTrueVSize() << endl;

    // Complex sesquilinear form and complex RHS
    SesquilinearForm sform(&fes);
    ComplexLinearForm rhs(&fes);

    // ComplexGridFunction to hold Dirichlet boundary (incident) values
    ComplexGridFunction x(&fes);
    x = 0.0;

    // Essential BCs: PEC on walls (attributes 1,2,5,6) and driven port 3
    Array<int> ess_bdr(mesh.bdr_attributes.Max()); ess_bdr = 0;
    // mark typical PEC faces (adjust to your mesh if different)
    if (ess_bdr.Size() >= 6)
    {
        ess_bdr[0] = 1; // attr 1
        ess_bdr[1] = 1; // attr 2
        ess_bdr[4] = 1; // attr 5
        ess_bdr[5] = 1; // attr 6
    }
    // drive port (attribute 3)
    if (ess_bdr.Size() >= 3) ess_bdr[2] = 1;

    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    // Project modal excitation on port 3: real TE10, imag = 0
    TE10Modal te10(a, 1.0);
    Array<int> port3(mesh.bdr_attributes.Max()); port3 = 0; if (port3.Size()>=3) port3[2]=1;
    x.real().ProjectBdrCoefficientTangent(te10, port3);
    // imag part zero (already zero)

    // Domain integrators: curl-curl (mu^{-1}) and mass (-k0^2 eps)
    ConstantCoefficient muinv(1.0 / mu_r);
    ConstantCoefficient masscoef( - (k0*k0) * eps_r );
    // Provide explicit zero-valued imaginary integrators so real/imag sparsity
    // patterns match (needed by ComplexUMFPackSolver).
    ConstantCoefficient zero_coeff(0.0);
    sform.AddDomainIntegrator(new CurlCurlIntegrator(muinv),
                             new CurlCurlIntegrator(zero_coeff));
    sform.AddDomainIntegrator(new VectorFEMassIntegrator(masscoef),
                             new VectorFEMassIntegrator(zero_coeff));

    // Modal termination on port 4 (attribute 4)
    Array<int> port4(mesh.bdr_attributes.Max()); port4 = 0; if (port4.Size()>=4) port4[3]=1;
    const double mu0 = 4.0*M_PI*1e-7;
    complex<double> kz = sqrt(complex<double>(k0*k0 - (M_PI/a)*(M_PI/a)));
    complex<double> Yc = kz / (omega * mu0); // modal admittance
    double Y_real = Yc.real();
    double Y_imag = Yc.imag();
    ModalMatrixCoefficient modalReal(3, te10, Y_real);
    ModalMatrixCoefficient modalImag(3, te10, Y_imag);
    sform.AddBoundaryIntegrator(new MixedVectorMassIntegrator(modalReal),
                                new MixedVectorMassIntegrator(modalImag),
                                port4);

    // Assemble and finalize (0 = integration order -> auto)
    sform.Assemble(0); sform.Finalize(0);
    rhs = 0.0; rhs.Assemble();

    // Form linear system and solve (GMRES on the returned operator)
    OperatorHandle Aop;
    Vector X, Bvec;
    sform.FormLinearSystem(ess_tdof_list, x, rhs, Aop, X, Bvec);

    cout << "Assembled operator: " << Aop.Type()
         << ", RHS norm = " << Bvec.Norml2() << endl;

#ifdef MFEM_USE_SUITESPARSE
    // Prefer direct complex UMFPACK when available
    ComplexSparseMatrix *csm = Aop.As<ComplexSparseMatrix>();
    if (csm)
    {
        cout << "Using ComplexUMFPack direct solver..." << endl;
        ComplexUMFPackSolver csolver(*csm);
        csolver.Mult(Bvec, X);
    }
    else
    {
        // Fallback to iterative solver if operator not assembled as ComplexSparseMatrix
        GMRESSolver solver;
        solver.SetOperator(*Aop.Ptr());
        solver.SetRelTol(1e-6);
        solver.SetAbsTol(1e-9);
        solver.SetMaxIter(2000);
        solver.SetKDim(min(100, (int)fes.GetTrueVSize()));
        solver.SetPrintLevel(0);
        solver.Mult(Bvec, X);
    }
#else
    GMRESSolver solver;
    solver.SetOperator(*Aop.Ptr());
    solver.SetRelTol(1e-6);
    solver.SetAbsTol(1e-9);
    solver.SetMaxIter(2000);
    solver.SetKDim(min(100, (int)fes.GetTrueVSize()));
    solver.SetPrintLevel(0);
    solver.Mult(Bvec, X);
#endif

    // Recover solution
    sform.RecoverFEMSolution(X, rhs, x);
    cout << "Solution recovered (real L2: " << x.real().Norml2()
         << ", imag L2: " << x.imag().Norml2() << ")" << endl;

    // Simple modal projection on a port: <emode, E> / <emode, emode>
    auto project_mode_on_port = [&](Array<int> &port_attr)->complex<real_t>
    {
        complex<real_t> num = 0.0;
        double den = 0.0;
        for (int be=0; be < fes.GetNBE(); ++be)
        {
            int attr = mesh.GetBdrAttribute(be);
            if (attr < 1 || attr > port_attr.Size() || port_attr[attr-1]==0) continue;
            const FiniteElement *bfe = fes.GetBE(be);
            ElementTransformation *Tr = fes.GetBdrElementTransformation(be);
            int qorder = max(2, 2*order+2);
            const IntegrationRule &ir = IntRules.Get(Tr->GetGeometryType(), qorder);
            for (int ipi=0; ipi < ir.GetNPoints(); ++ipi)
            {
                const IntegrationPoint &ip = ir.IntPoint(ipi);
                Tr->SetIntPoint(&ip);
                Vector em(3); te10.Eval(em, *Tr, ip);
                Vector Er(3), Ei(3);
                x.real().GetVectorValue(*Tr, ip, Er);
                x.imag().GetVectorValue(*Tr, ip, Ei);
                double dot_r = em * Er;
                double dot_i = em * Ei;
                complex<real_t> dot(dot_r, dot_i);
                double w = ip.weight * Tr->Weight();
                num += dot * w;
                den += (em * em) * w;
            }
        }
        if (den == 0.0) return complex<real_t>(0.0);
        return num / den;
    };

    Array<int> p3 = port3; Array<int> p4 = port4;
    complex<real_t> a_tot3 = project_mode_on_port(p3);
    complex<real_t> a_tot4 = project_mode_on_port(p4);
    complex<real_t> a_inc(1.0, 0.0);
    complex<real_t> a_ref = a_tot3 - a_inc;
    complex<real_t> S11 = (a_inc != complex<real_t>(0.0)) ? (a_ref / a_inc) : complex<real_t>(0.0);
    complex<real_t> S21 = (a_inc != complex<real_t>(0.0)) ? (a_tot4 / a_inc) : complex<real_t>(0.0);

    cout << "S11 = " << S11 << ", S21 = " << S21 << ", |S21| = " << abs(S21) << endl;

    // Save fields for visualization (separate real/imag)
    ParaViewDataCollection paraview(output_dir, &mesh);
    paraview.RegisterField("E_real", &x.real());
    paraview.RegisterField("E_imag", &x.imag());
    paraview.SetDataFormat(mfem::VTKFormat::ASCII);
    paraview.Save();

    // Write S-params
    ofstream of((string(output_dir)+"/s_parameters.txt").c_str());
    of << "S11 " << S11.real() << " " << S11.imag() << "\n";
    of << "S21 " << S21.real() << " " << S21.imag() << "\n";
    of.close();

    return 0;
}
