#include "mfem.hpp"
#include "hpcfem/solvers/solver_hypre_amg.hpp"
#include <mpi.h>
#include <fstream>
#include <iostream>
#include <map>
#include <cmath>
#include <chrono>

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
    // Read simple command-line args before MPI_Init to ensure they're
    // available on all MPI implementations (some may modify argc/argv).
    const char *mesh_file = "testdata/testmesh_package.mesh";
    const char *output_dir = "results/ex6_multiple_materials";
    int order = 2;

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

    // Initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm comm = MPI_COMM_WORLD;
    int myid = 0;
    int numprocs = 1;
    MPI_Comm_rank(comm, &myid);
    MPI_Comm_size(comm, &numprocs);
    // ============================================================
    // STEP 1: PROBLEM SETUP
    // ============================================================

    if (myid == 0)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Multi-Material Thermal Analysis" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }

    
    // STEP 2: LOAD MESH
    // ============================================================

    // We'll perform preprocessing (mesh I/O, FE space setup and serial assembly)
    // on rank 0 only. The solver will run in parallel across all ranks.

    // For Option B we perform mesh reading and preprocessing in parallel.
    // Create the serial Mesh object on each rank and then construct the
    // ParMesh by partitioning it. We keep I/O/prints on rank 0 only.
    Mesh mesh(mesh_file, 1, 1);
    if (myid == 0)
    {
        std::cout << "\nMesh information:" << std::endl;
        std::cout << "  Dimension:     " << mesh.Dimension() << std::endl;
        std::cout << "  Elements:      " << mesh.GetNE() << std::endl;
        std::cout << "  Boundary elem: " << mesh.GetNBE() << std::endl;

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
    }

    ParMesh *pmesh = new ParMesh(MPI_COMM_WORLD, mesh);
    int dim = pmesh->Dimension();

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

    if (myid == 0)
    {
        std::cout << "\nMaterial properties (converted to mm units):" << std::endl;
        std::cout << "  k_FR4:      " << k_FR4 << " W/mm/K" << std::endl;
        std::cout << "  k_Plastic:  " << k_Plastic << " W/mm/K" << std::endl;
        std::cout << "  k_Silicon:  " << k_Silicon << " W/mm/K" << std::endl;
        std::cout << "  k_Aluminum: " << k_Aluminum << " W/mm/K" << std::endl;
        std::cout << "  Q_Silicon:  " << Q_Silicon << " W/mm³" << std::endl;
        std::cout << "  h_conv:     " << h_conv << " W/mm²/K" << std::endl;
        std::cout << "  T_amb:      " << T_amb << " K" << std::endl;
        std::cout << "  T_BC:       " << T_dirichlet << " K" << std::endl;
    }

    // ============================================================
    // STEP 4: FINITE ELEMENT SPACE
    // ============================================================

    H1_FECollection fec(order, dim);
    ParFiniteElementSpace fes(pmesh, &fec);

    if (myid == 0)
    {
        std::cout << "\nFinite element space:" << std::endl;
        std::cout << "  Number of DOFs: " << fes.GetTrueVSize() << std::endl;
    }

    // ============================================================
    // STEP 5: DEFINE COEFFICIENTS FOR MATERIALS
    // ============================================================

    // Thermal conductivity coefficient (piecewise constant by domain)
    PiecewiseConstantCoefficient k_coef(pmesh, k_Aluminum); // default is Aluminum
    
    // Map domain attributes to materials
    k_coef.SetMaterialProperty(1, k_FR4);      // Domain 1 (FR4)
    k_coef.SetMaterialProperty(2, k_Plastic); // Domain 2 (Plastic)
    k_coef.SetMaterialProperty(11, k_Silicon); // Domain 11 (Silicon)

    if (myid == 0)
    {
        std::cout << "\nMaterial assignment (domain attribute → material):" << std::endl;
        std::cout << "  Attribute 1:  FR4 (k = " << k_FR4 << " W/mm/K)" << std::endl;
        std::cout << "  Attribute 2: Plastic (k = " << k_Plastic << " W/mm/K)" << std::endl;
        std::cout << "  Attribute 11: Silicon (k = " << k_Silicon << " W/mm/K)" << std::endl;
        std::cout << "  All others:   Aluminum (k = " << k_Aluminum << " W/mm/K)" << std::endl;
    }

    // Volumetric heat source (only in Silicon domain)
    PiecewiseSourceCoefficient Q_coef(pmesh, 0.0); // default is no source
    Q_coef.SetSourceInDomain(11, Q_Silicon);       // Heat source in Silicon (attribute 11)

    if (myid == 0)
    {
        std::cout << "\nHeat source:" << std::endl;
        std::cout << "  Domain 11 (Silicon): Q = " << Q_Silicon << " W/mm³" << std::endl;
    }

    // ============================================================
    // STEP 6: BOUNDARY CONDITIONS
    // ============================================================

    // Identify exterior boundaries (parallel-safe): mark a boundary attribute as
    // Robin (exterior) if any local boundary element with that attribute is an
    // exterior face on any rank. We build a local marker and do an MPI_Allreduce
    // with MPI_MAX so every rank has the global marker.
    int local_max_bdr_attr = pmesh->bdr_attributes.Max();
    int global_max_bdr_attr = 0;
    MPI_Allreduce(&local_max_bdr_attr, &global_max_bdr_attr, 1, MPI_INT, MPI_MAX, comm);

    Array<int> robin_bdr(global_max_bdr_attr);
    robin_bdr = 0;

    for (int be = 0; be < pmesh->GetNBE(); be++)
    {
        int bdr_attr = pmesh->GetBdrAttribute(be);
        int face_idx, orientation;
        pmesh->GetBdrElementFace(be, &face_idx, &orientation);
        if (!pmesh->FaceIsInterior(face_idx))
        {
            if (bdr_attr > 0 && bdr_attr <= robin_bdr.Size())
            {
                robin_bdr[bdr_attr - 1] = 1;
            }
        }
    }

    // Combine across ranks so every rank has the full list
    if (global_max_bdr_attr > 0)
    {
        std::vector<int> local_vec(global_max_bdr_attr, 0), global_vec(global_max_bdr_attr, 0);
        for (int i = 0; i < robin_bdr.Size(); ++i) local_vec[i] = robin_bdr[i];
        MPI_Allreduce(local_vec.data(), global_vec.data(), global_max_bdr_attr, MPI_INT, MPI_MAX, comm);
        for (int i = 0; i < global_max_bdr_attr; ++i) robin_bdr[i] = global_vec[i];
    }

    if (myid == 0)
    {
        std::cout << "  Exterior boundary attributes (robin): ";
        bool first = true;
        for (int i = 0; i < robin_bdr.Size(); ++i)
        {
            if (robin_bdr[i])
            {
                if (!first) std::cout << ", ";
                std::cout << (i + 1);
                first = false;
            }
        }
        std::cout << std::endl;
    }

    // Dirichlet BC: boundary 4 (1-based in mesh file, but 0-based in MFEM)
    Array<int> ess_bdr(global_max_bdr_attr);
    ess_bdr = 0;
    if (global_max_bdr_attr >= 4)
    {
        ess_bdr[4 - 1] = 1; // Boundary 4 (0-based index: 3)
    }

    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    if (myid == 0)
    {
        std::cout << "\nBoundary conditions:" << std::endl;
        std::cout << "  Dirichlet BC: boundary 4, T = " << T_dirichlet << " K" << std::endl;
        std::cout << "  Essential DOFs: " << ess_tdof_list.Size() << std::endl;
    }

    // Remove Dirichlet boundary from Robin markers
    if (global_max_bdr_attr >= 4)
    {
        robin_bdr[4 - 1] = 0; // Remove boundary 4 (Dirichlet)
    }

    if (myid == 0)
    {
        std::cout << "  Robin BC (convection): applied only to exterior boundaries" << std::endl;
        std::cout << "    h = " << h_conv << " W/mm²/K, T_amb = " << T_amb << " K" << std::endl;
        std::cout << "    (Interior boundaries excluded)" << std::endl;
    }

    // ============================================================
    // STEP 7: ASSEMBLE THE LINEAR SYSTEM
    // ============================================================

    if (myid == 0)
    {
        std::cout << "\nAssembling linear system..." << std::endl;
    }

    std::chrono::high_resolution_clock::time_point solving_start = std::chrono::high_resolution_clock::now();

    // Bilinear form: a(u,v) = ∫ k ∇u·∇v dx + ∫ h u·v ds
    ParBilinearForm a(&fes);
    a.AddDomainIntegrator(new DiffusionIntegrator(k_coef));

    // Add Robin BC (convection term): ∫ h u v ds
    ConstantCoefficient h_coef(h_conv);
    a.AddBoundaryIntegrator(new MassIntegrator(h_coef), robin_bdr);

    a.Assemble();
    a.Finalize();

    // Linear form: L(v) = ∫ Q v dx + ∫ h T_amb v ds
    ParLinearForm b(&fes);
    b.AddDomainIntegrator(new DomainLFIntegrator(Q_coef));

    // Robin BC RHS: ∫ h T_amb v ds
    ConstantCoefficient hT_amb(h_conv * T_amb);
    b.AddBoundaryIntegrator(new BoundaryLFIntegrator(hT_amb), robin_bdr);

    b.Assemble();

    if (myid == 0)
    {
        std::cout << "  System assembled." << std::endl;
    }

    std::chrono::high_resolution_clock::time_point assembly_end = std::chrono::high_resolution_clock::now();

    // ============================================================
    // STEP 8: APPLY DIRICHLET BOUNDARY CONDITIONS
    // ============================================================

    ParGridFunction T(&fes);
    T = T_amb; // Initial guess

    // Project Dirichlet boundary values
    ConstantCoefficient T_D(T_dirichlet);
    T.ProjectBdrCoefficient(T_D, ess_bdr);

    // ============================================================
    // STEP 9: SOLVE THE LINEAR SYSTEM
    // ============================================================

    if (myid == 0)
    {
        std::cout << "\nSolving linear system..." << std::endl;
    }

    HypreParMatrix A;
    Vector B, X;
    a.FormLinearSystem(ess_tdof_list, T, b, A, X, B);

    if (myid == 0)
    {
        // Print global system size (number of true DOFs)
        HYPRE_BigInt global_rows = fes.GlobalTrueVSize();
        std::cout << "  System global true DOFs: " << global_rows << std::endl;
    }

    // Solve with Hypre AMG wrapper
    hpcfem::HypreAmgSolver amg_solver(1.0e-9, 1000, 0);
    amg_solver.solve(A, B, X);

    if (myid == 0)
    {
        std::cout << "  Hypre AMG iterations: " << amg_solver.getNumIterations() << std::endl;
        std::cout << "  Hypre AMG final norm: " << amg_solver.getFinalNorm() << std::endl;
    }

    // Recover the solution
    a.RecoverFEMSolution(X, b, T);

   // ============================================================
    // STEP 10: CHECK SOLUTION AND OUTPUT RESULTS
    // ============================================================

    // Compute global min/max/avg via explicit MPI reductions so the values
    // reflect the final (post-projection) ParGridFunction across all ranks.
    double local_min = std::numeric_limits<double>::infinity();
    double local_max = -std::numeric_limits<double>::infinity();
    double local_sum = 0.0;
    int local_count = 0;
    for (int i = 0; i < T.Size(); ++i)
    {
        double v = T(i);
        if (v < local_min) local_min = v;
        if (v > local_max) local_max = v;
        local_sum += v;
        ++local_count;
    }
    double T_min = 0.0, T_max = 0.0, T_avg = 0.0;
    double global_sum = 0.0;
    int global_count = 0;
    MPI_Allreduce(&local_min, &T_min, 1, MPI_DOUBLE, MPI_MIN, comm);
    MPI_Allreduce(&local_max, &T_max, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, comm);
    MPI_Allreduce(&local_count, &global_count, 1, MPI_INT, MPI_SUM, comm);
    if (global_count > 0) T_avg = global_sum / static_cast<double>(global_count);

    std::chrono::high_resolution_clock::time_point solving_end = std::chrono::high_resolution_clock::now();

    if (myid == 0)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "SOLUTION SUMMARY" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Temperature range:" << std::endl;
        std::cout << "  Min: " << T_min << " K" << std::endl;
        std::cout << "  Max: " << T_max << " K" << std::endl;
        std::cout << "  Avg: " << T_avg << " K" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Total time: " << std::chrono::duration_cast<std::chrono::milliseconds>(solving_end - solving_start).count() << " ms" << std::endl;
        std::cout << "Assembly time: " << std::chrono::duration_cast<std::chrono::milliseconds>(assembly_end - solving_start).count() << " ms" << std::endl;
        std::cout << "Solving time: " << std::chrono::duration_cast<std::chrono::milliseconds>(solving_end - assembly_end).count() << " ms" << std::endl;
    }

    // Sanity check
    const double T_expected_min = 303.0;
    const double T_expected_max = 325.0;
    
    if (T_min < T_expected_min || T_max > T_expected_max)
    {
        if (myid == 0)
        {
            std::cout << "\nWARNING: Temperature out of expected range!" << std::endl;
            std::cout << "   Expected range: [" << T_expected_min << ", " 
                      << T_expected_max << "] K" << std::endl;
            std::cout << "   Actual range:   [" << T_min << ", " << T_max << "] K" << std::endl;
        }
    }
    else
    {
        if (myid == 0)
        {
            std::cout << "\n✓ Temperature range is physically reasonable." << std::endl;
        }
    }

    // ============================================================
    // STEP 11: SAVE SOLUTION FOR VISUALIZATION
    // ============================================================

    if (myid == 0)
    {
        std::cout << "\nSaving solution to " << output_dir << " ..." << std::endl;
    }

    ParaViewDataCollection paraview_dc(output_dir, pmesh);
    paraview_dc.SetLevelsOfDetail(order);
    paraview_dc.RegisterField("Temperature", &T);
    paraview_dc.SetTime(0.0);
    paraview_dc.SetCycle(0);
    paraview_dc.Save();

    if (myid == 0)
    {
        std::cout << "Solution saved successfully." << std::endl;
        std::cout << "\nTo visualize:" << std::endl;
        std::cout << "  paraview " << output_dir << "/ex6_multiple_materials.pvd" << std::endl;

        std::cout << "\n========================================" << std::endl;
        std::cout << "Simulation completed successfully!" << std::endl;
        std::cout << "========================================\n" << std::endl;
    }
    // Finalize MPI
    MPI_Finalize();

    return 0;
}
