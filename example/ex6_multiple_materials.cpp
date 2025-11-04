#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <map>
#include <cmath>

/**
 * Example 6: Multi-Material Thermal Analysis
 *
 * A stationary (steady-state) thermal problem with multiple materials.
 * Different materials have different thermal conductivities.
 *
 * PROBLEM SETUP:
 * - Domain 1 :   FR4 (k = 0.3 W/m/K)
 * - Domain 2 :  Plastic (k = 0.2 W/m/K)
 * - Domain 11 : Silicon (k = 150 W/m/K)
 * - All other domains:        Aluminum (k = 205 W/m/K)
 *
 * BOUNDARY CONDITIONS:
 * - Volumetric heat source: Q = 2e8 W/m³ in Domain 11 (Silicon)
 *   (scaled to 2e2 W/mm³ due to mm units)
 * - Robin BC (convection): h = 50 W/m²/K, T_amb = 303.15 K on all exterior boundaries
 *   (scaled to 0.05 W/mm²/K due to mm units)
 * - Dirichlet BC: T = 323.15 K on boundary 4
 *
 * GOVERNING EQUATION (steady-state):
 *   -∇·(k(x)∇T) = Q(x)  in Ω
 *   T = T_D            on Γ_D (Dirichlet boundary)
 *   -k∇T·n = h(T - T_amb) on Γ_N (Robin boundary - convection)
 *
 * UNITS:
 * - Length: mm (mesh is in mm)
 * - Temperature: K
 * - Power density: W/mm³
 * - Conductivity: W/mm/K
 * - Convection: W/mm²/K
 */

using namespace mfem;

// Piecewise constant coefficient for thermal conductivity based on element attributes
class PiecewiseConstantCoefficient : public Coefficient
{
private:
    std::map<int, double> attr_to_value_;
    double default_value_;
    Mesh *mesh_;

public:
    PiecewiseConstantCoefficient(Mesh *mesh, double default_val = 1.0)
        : mesh_(mesh), default_value_(default_val) {}

    void SetMaterialProperty(int attr, double value)
    {
        attr_to_value_[attr] = value;
    }

    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        int attr = T.Attribute;
        auto it = attr_to_value_.find(attr);
        if (it != attr_to_value_.end())
        {
            return it->second;
        }
        return default_value_;
    }
};

// Piecewise constant source term based on element attributes
class PiecewiseSourceCoefficient : public Coefficient
{
private:
    std::map<int, double> attr_to_value_;
    double default_value_;
    Mesh *mesh_;

public:
    PiecewiseSourceCoefficient(Mesh *mesh, double default_val = 0.0)
        : mesh_(mesh), default_value_(default_val) {}

    void SetSourceInDomain(int attr, double value)
    {
        attr_to_value_[attr] = value;
    }

    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        int attr = T.Attribute;
        auto it = attr_to_value_.find(attr);
        if (it != attr_to_value_.end())
        {
            return it->second;
        }
        return default_value_;
    }
};

int main(int argc, char *argv[])
{
    // ============================================================
    // STEP 1: PROBLEM SETUP
    // ============================================================

    std::cout << "\n========================================" << std::endl;
    std::cout << "Multi-Material Thermal Analysis" << std::endl;
    std::cout << "========================================\n" << std::endl;

    const char *mesh_file = "testdata/testmesh_package.mesh";
    const char *output_dir = "results/ex6_multiple_materials";
    int order = 1;

    if (argc > 1)
    {
        mesh_file = argv[1];
    }
    if (argc > 2)
    {
        output_dir = argv[2];
    }
    if (argc > 3)
    {
        order = atoi(argv[3]);
    }

    std::cout << "Mesh file:     " << mesh_file << std::endl;
    std::cout << "Output dir:    " << output_dir << std::endl;
    std::cout << "FE order:      " << order << std::endl;

    // ============================================================
    // STEP 2: LOAD MESH
    // ============================================================

    Mesh mesh(mesh_file);
    int dim = mesh.Dimension();

    std::cout << "\nMesh information:" << std::endl;
    std::cout << "  Dimension:     " << dim << std::endl;
    std::cout << "  Elements:      " << mesh.GetNE() << std::endl;
    std::cout << "  Boundary elem: " << mesh.GetNBE() << std::endl;

    // Print element attributes (domain IDs)
    std::cout << "\nElement attributes (domain IDs):" << std::endl;
    for (int i = 0; i < mesh.attributes.Size(); i++)
    {
        int attr = mesh.attributes[i];
        int count = 0;
        for (int e = 0; e < mesh.GetNE(); e++)
        {
            if (mesh.GetAttribute(e) == attr)
                count++;
        }
        std::cout << "  Attribute " << attr << ": " << count << " elements" << std::endl;
    }

    // ============================================================
    // STEP 3: MATERIAL PROPERTIES (in mm units)
    // ============================================================
    // NOTE: mesh is in mm, so we need to scale units accordingly
    // Original: W/m/K → W/mm/K (divide by 1000)
    // Original: W/m³ → W/mm³ (divide by 1e9)
    // Original: W/m²/K → W/mm²/K (divide by 1000)

    const double SCALE_LENGTH = 1000.0; // m to mm

    // Thermal conductivities [W/mm/K]
    const double k_FR4 = 0.3 / SCALE_LENGTH;         // FR4: 0.3 W/m/K
    const double k_Plastic = 0.2 / SCALE_LENGTH;     // Plastic: 0.2 W/m/K
    const double k_Silicon = 150.0 / SCALE_LENGTH;   // Silicon: 150 W/m/K
    const double k_Aluminum = 205.0 / SCALE_LENGTH;  // Aluminum: 205 W/m/K

    // Heat source [W/mm³]
    const double Q_Silicon = 2e8 / (SCALE_LENGTH * SCALE_LENGTH * SCALE_LENGTH); // 2e8 W/m³

    // Convection coefficient [W/mm²/K]
    const double h_conv = 50.0 / (SCALE_LENGTH * SCALE_LENGTH); // 50 W/m²/K

    // Temperatures [K]
    const double T_amb = 303.15;      // Ambient temperature
    const double T_dirichlet = 323.15; // Dirichlet BC temperature

    std::cout << "\nMaterial properties (converted to mm units):" << std::endl;
    std::cout << "  k_FR4:      " << k_FR4 << " W/mm/K" << std::endl;
    std::cout << "  k_Plastic:  " << k_Plastic << " W/mm/K" << std::endl;
    std::cout << "  k_Silicon:  " << k_Silicon << " W/mm/K" << std::endl;
    std::cout << "  k_Aluminum: " << k_Aluminum << " W/mm/K" << std::endl;
    std::cout << "  Q_Silicon:  " << Q_Silicon << " W/mm³" << std::endl;
    std::cout << "  h_conv:     " << h_conv << " W/mm²/K" << std::endl;
    std::cout << "  T_amb:      " << T_amb << " K" << std::endl;
    std::cout << "  T_BC:       " << T_dirichlet << " K" << std::endl;

    // ============================================================
    // STEP 4: FINITE ELEMENT SPACE
    // ============================================================

    H1_FECollection fec(order, dim);
    FiniteElementSpace fes(&mesh, &fec);

    std::cout << "\nFinite element space:" << std::endl;
    std::cout << "  Number of DOFs: " << fes.GetTrueVSize() << std::endl;

    // ============================================================
    // STEP 5: DEFINE COEFFICIENTS FOR MATERIALS
    // ============================================================

    // Thermal conductivity coefficient (piecewise constant by domain)
    PiecewiseConstantCoefficient k_coef(&mesh, k_Aluminum); // default is Aluminum
    
    // Map domain attributes to materials
    k_coef.SetMaterialProperty(1, k_FR4);      // Domain 1 (FR4)
    k_coef.SetMaterialProperty(2, k_Plastic); // Domain 2 (Plastic)
    k_coef.SetMaterialProperty(11, k_Silicon); // Domain 11 (Silicon)

    std::cout << "\nMaterial assignment (domain attribute → material):" << std::endl;
    std::cout << "  Attribute 1:  FR4 (k = " << k_FR4 << " W/mm/K)" << std::endl;
    std::cout << "  Attribute 2: Plastic (k = " << k_Plastic << " W/mm/K)" << std::endl;
    std::cout << "  Attribute 11: Silicon (k = " << k_Silicon << " W/mm/K)" << std::endl;
    std::cout << "  All others:   Aluminum (k = " << k_Aluminum << " W/mm/K)" << std::endl;

    // Volumetric heat source (only in Silicon domain)
    PiecewiseSourceCoefficient Q_coef(&mesh, 0.0); // default is no source
    Q_coef.SetSourceInDomain(11, Q_Silicon);       // Heat source in Silicon (attribute 11)

    std::cout << "\nHeat source:" << std::endl;
    std::cout << "  Domain 11 (Silicon): Q = " << Q_Silicon << " W/mm³" << std::endl;

    // ============================================================
    // STEP 6: BOUNDARY CONDITIONS
    // ============================================================

    // Identify exterior boundaries using MFEM's FaceIsInterior() method
    // Interior boundaries (material interfaces) should NOT have convection BC
    std::cout << "\nIdentifying exterior vs interior boundaries..." << std::endl;
    
    Array<int> exterior_bdr_attrs;
    Array<int> interior_bdr_attrs;
    
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
        int bdr_attr = mesh.GetBdrAttribute(be);
        int face_idx, orientation;
        mesh.GetBdrElementFace(be, &face_idx, &orientation);
        
        // Check if this face is interior (shared between elements)
        if (mesh.FaceIsInterior(face_idx))
        {
            // This is an interior boundary (material interface)
            if (interior_bdr_attrs.Find(bdr_attr) == -1)
            {
                interior_bdr_attrs.Append(bdr_attr);
            }
        }
        else
        {
            // This is an exterior boundary
            if (exterior_bdr_attrs.Find(bdr_attr) == -1)
            {
                exterior_bdr_attrs.Append(bdr_attr);
            }
        }
    }
    
    // Sort for cleaner output
    exterior_bdr_attrs.Sort();
    interior_bdr_attrs.Sort();
    
    std::cout << "  Exterior boundaries: ";
    for (int i = 0; i < exterior_bdr_attrs.Size(); i++)
    {
        std::cout << exterior_bdr_attrs[i];
        if (i < exterior_bdr_attrs.Size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;
    
    std::cout << "  Interior boundaries: ";
    for (int i = 0; i < interior_bdr_attrs.Size(); i++)
    {
        std::cout << interior_bdr_attrs[i];
        if (i < interior_bdr_attrs.Size() - 1) std::cout << ", ";
    }
    std::cout << std::endl;

    // Dirichlet BC: boundary 4 (1-based in mesh file, but 0-based in MFEM)
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 0;
    if (mesh.bdr_attributes.Max() >= 4)
    {
        ess_bdr[4 - 1] = 1; // Boundary 4 (0-based index: 3)
    }

    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    std::cout << "\nBoundary conditions:" << std::endl;
    std::cout << "  Dirichlet BC: boundary 4, T = " << T_dirichlet << " K" << std::endl;
    std::cout << "  Essential DOFs: " << ess_tdof_list.Size() << std::endl;

    // Robin BC (convection): only on exterior boundaries, excluding Dirichlet BC
    Array<int> robin_bdr(mesh.bdr_attributes.Max());
    robin_bdr = 0; // Start with no boundaries
    
    // Apply Robin BC only to exterior boundaries
    for (int i = 0; i < exterior_bdr_attrs.Size(); i++)
    {
        int attr = exterior_bdr_attrs[i];
        if (attr <= mesh.bdr_attributes.Max())
        {
            robin_bdr[attr - 1] = 1; // 0-based indexing
        }
    }
    
    // Remove boundary 4 (Dirichlet BC)
    if (mesh.bdr_attributes.Max() >= 4)
    {
        robin_bdr[4 - 1] = 0; // Remove boundary 4 (Dirichlet)
    }

    std::cout << "  Robin BC (convection): applied only to exterior boundaries" << std::endl;
    std::cout << "    h = " << h_conv << " W/mm²/K, T_amb = " << T_amb << " K" << std::endl;
    std::cout << "    (Interior boundaries excluded)" << std::endl;

    // ============================================================
    // STEP 7: ASSEMBLE THE LINEAR SYSTEM
    // ============================================================

    std::cout << "\nAssembling linear system..." << std::endl;

    // Bilinear form: a(u,v) = ∫ k ∇u·∇v dx + ∫ h u·v ds
    BilinearForm a(&fes);
    a.AddDomainIntegrator(new DiffusionIntegrator(k_coef));
    
    // Add Robin BC (convection term): ∫ h u v ds
    ConstantCoefficient h_coef(h_conv);
    a.AddBoundaryIntegrator(new MassIntegrator(h_coef), robin_bdr);
    
    a.Assemble();
    a.Finalize();

    // Linear form: L(v) = ∫ Q v dx + ∫ h T_amb v ds
    LinearForm b(&fes);
    b.AddDomainIntegrator(new DomainLFIntegrator(Q_coef));
    
    // Robin BC RHS: ∫ h T_amb v ds
    ConstantCoefficient hT_amb(h_conv * T_amb);
    b.AddBoundaryIntegrator(new BoundaryLFIntegrator(hT_amb), robin_bdr);
    
    b.Assemble();

    std::cout << "  System assembled." << std::endl;

    // ============================================================
    // STEP 8: APPLY DIRICHLET BOUNDARY CONDITIONS
    // ============================================================

    GridFunction T(&fes);
    T = T_amb; // Initial guess

    // Project Dirichlet boundary values
    ConstantCoefficient T_D(T_dirichlet);
    T.ProjectBdrCoefficient(T_D, ess_bdr);

    // ============================================================
    // STEP 9: SOLVE THE LINEAR SYSTEM
    // ============================================================

    std::cout << "\nSolving linear system..." << std::endl;

    SparseMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, T, b, A, X, B);

    std::cout << "  System size: " << A.Height() << " x " << A.Width() << std::endl;

    // Use CG solver with GS preconditioner
    GSSmoother M(A);
    CGSolver solver;
    solver.SetRelTol(1e-9);
    solver.SetMaxIter(10000);
    solver.SetPrintLevel(1);
    solver.SetPreconditioner(M);
    solver.SetOperator(A);

    solver.Mult(B, X);

    if (solver.GetConverged())
    {
        std::cout << "  PCG solver converged in " << solver.GetNumIterations()
                  << " iterations." << std::endl;
    }
    else
    {
        std::cout << "  WARNING: PCG solver did NOT converge!" << std::endl;
    }

    // Recover the solution
    a.RecoverFEMSolution(X, b, T);

    // ============================================================
    // STEP 10: CHECK SOLUTION AND OUTPUT RESULTS
    // ============================================================

    double T_min = T.Min();
    double T_max = T.Max();
    double T_avg = T.Sum() / T.Size();

    std::cout << "\n========================================" << std::endl;
    std::cout << "SOLUTION SUMMARY" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Temperature range:" << std::endl;
    std::cout << "  Min: " << T_min << " K" << std::endl;
    std::cout << "  Max: " << T_max << " K" << std::endl;
    std::cout << "  Avg: " << T_avg << " K" << std::endl;
    std::cout << "========================================" << std::endl;

    // Sanity check
    const double T_expected_min = 303.0;
    const double T_expected_max = 325.0;
    
    if (T_min < T_expected_min || T_max > T_expected_max)
    {
        std::cout << "\nWARNING: Temperature out of expected range!" << std::endl;
        std::cout << "   Expected range: [" << T_expected_min << ", " 
                  << T_expected_max << "] K" << std::endl;
        std::cout << "   Actual range:   [" << T_min << ", " << T_max << "] K" << std::endl;
    }
    else
    {
        std::cout << "\n✓ Temperature range is physically reasonable." << std::endl;
    }

    // ============================================================
    // STEP 11: SAVE SOLUTION FOR VISUALIZATION
    // ============================================================

    std::cout << "\nSaving solution to " << output_dir << " ..." << std::endl;

    ParaViewDataCollection paraview_dc(output_dir, &mesh);
    paraview_dc.SetLevelsOfDetail(order);
    paraview_dc.RegisterField("Temperature", &T);
    paraview_dc.SetTime(0.0);
    paraview_dc.SetCycle(0);
    paraview_dc.Save();

    std::cout << "Solution saved successfully." << std::endl;
    std::cout << "\nTo visualize:" << std::endl;
    std::cout << "  paraview " << output_dir << "/ex6_multiple_materials.pvd" << std::endl;

    std::cout << "\n========================================" << std::endl;
    std::cout << "Simulation completed successfully!" << std::endl;
    std::cout << "========================================\n" << std::endl;

    return 0;
}
