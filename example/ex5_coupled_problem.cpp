#include "mfem.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

/**
 * Example: Electric-Thermal Coupling (Joule Heating)
 *
 * This demonstrates a fully coupled multiphysics problem where:
 * 1. Electric field generates Joule heating: Q = σ(T) |∇φ|²
 * 2. Temperature affects electrical conductivity: σ(T) = σ₀/(1 + α(T - T₀))
 *
 * GOVERNING EQUATIONS:
 * -----------------
 * Electric potential (elliptic):
 *   -∇·(σ(T)∇φ) = 0  in Ω
 *   φ = φ_prescribed on Γ_D (Dirichlet boundaries)
 *
 * Temperature (parabolic):
 *   ρc ∂T/∂t - ∇·(k∇T) = Q(φ,T)  in Ω
 *   T = T_prescribed on Γ_D
 *   -k∇T·n = h(T - T_amb) on Γ_N (Robin BC - convection)
 *
 * COUPLING PROCEDURE:
 * ------------------
 * At each time step n→n+1:
 * 1. Given T^n, compute σ(T^n)
 * 2. Solve electric problem: -∇·(σ(T^n)∇φ^{n+1}) = 0
 * 3. Compute Joule heating: Q^{n+1} = σ(T^n)|∇φ^{n+1}|²
 * 4. Solve thermal problem: M(T^{n+1} - T^n)/dt + K(T^n)T^{n+1} = b + Q^{n+1}
 * 5. Check convergence: if |T^{n+1} - T^n| < tol, continue; else iterate
 *
 * KEY MFEM OPERATIONS DEMONSTRATED:
 * ---------------------------------
 * - Multiple FiniteElementSpaces (one for φ, one for T)
 * - GridFunction evaluation and coefficient updates
 * - BilinearForm assembly with temperature-dependent coefficients
 * - LinearForm assembly with computed source terms
 * - Time stepping with operator splitting
 * - Picard iteration for nonlinear coupling
 */

using namespace mfem;

// Temperature-dependent electrical conductivity: σ(T) = σ₀/(1 + α(T - T₀))
class TemperatureDependentConductivity : public Coefficient
{
private:
    const GridFunction &temp_;
    double sigma_0_; // Base conductivity [S/m]
    double alpha_;   // Temperature coefficient [1/K]
    double T_ref_;   // Reference temperature [K]

public:
    TemperatureDependentConductivity(const GridFunction &T,
                                     double sigma0 = 1e7,
                                     double alpha = 0.004,
                                     double Tref = 293.15)
        : temp_(T), sigma_0_(sigma0), alpha_(alpha), T_ref_(Tref) {}

    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        double temp_val = temp_.GetValue(T, ip);
        return sigma_0_ / (1.0 + alpha_ * (temp_val - T_ref_));
    }
};

// Joule heating source: Q = σ(T)|∇φ|²
class JouleHeatingCoefficient : public Coefficient
{
private:
    const GridFunction &potential_;
    const GridFunction &temp_;
    double sigma_0_;
    double alpha_;
    double T_ref_;
    mutable Vector grad_phi_;

public:
    JouleHeatingCoefficient(const GridFunction &phi, const GridFunction &T,
                            double sigma0 = 1e7,
                            double alpha = 0.004,
                            double Tref = 293.15)
        : potential_(phi), temp_(T),
          sigma_0_(sigma0), alpha_(alpha), T_ref_(Tref) {}

    virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip)
    {
        // Get temperature
        double temp_val = temp_.GetValue(T, ip);

        // Compute conductivity at current temperature
        double sigma = sigma_0_ / (1.0 + alpha_ * (temp_val - T_ref_));

        // Get gradient of potential: ∇φ
        potential_.GetGradient(T, grad_phi_);

        // Compute |∇φ|²
        double grad_phi_sq = grad_phi_ * grad_phi_;

        // Q = σ|∇φ|²
        return sigma * grad_phi_sq;
    }
};

int main(int argc, char *argv[])
{
    // ============================================================
    // STEP 1: PROBLEM SETUP
    // ============================================================

    std::cout << "\n========================================" << std::endl;
    std::cout << "Electric-Thermal Coupling Example" << std::endl;
    std::cout << "========================================\n"
              << std::endl;

    // Mesh parameters
    int nx = 20, ny = 10, nz = 10;
    int order = 1;

    // Physical parameters
    const double sigma_0 = 5.96e7;     // Copper conductivity at 20°C [S/m]
    const double alpha_elec = 0.00393; // Copper temp coefficient [1/K]
    const double T_ref = 293.15;       // Reference temperature [K]

    const double k_thermal = 400.0; // Thermal conductivity [W/m/K]
    const double rho = 8960.0;      // Density [kg/m³]
    const double cp = 385.0;        // Specific heat [J/kg/K]
    const double h_conv = 100.0;     // Convection coefficient [W/m²/K]
    const double T_amb = 293.15;    // Ambient temperature [K]

    // Applied voltage
    const double V_applied = 0.1; // [V]

    // Time stepping
    const double t_final = 100.0; // [s]
    const double dt = 1.0;        // [s]

    // Coupling iteration parameters
    const int max_coupling_iter = 20;
    const double coupling_tol = 1e-3;

    std::cout << "Physical parameters:" << std::endl;
    std::cout << "  Electrical conductivity: " << sigma_0 << " S/m" << std::endl;
    std::cout << "  Thermal conductivity:    " << k_thermal << " W/m/K" << std::endl;
    std::cout << "  Applied voltage:         " << V_applied << " V" << std::endl;
    std::cout << "  Time step:               " << dt << " s" << std::endl;
    std::cout << "========================================\n"
              << std::endl;

    // ============================================================
    // STEP 2: MESH GENERATION
    // ============================================================

    // Create rectangular bar mesh (conductor)
    Mesh mesh = Mesh::MakeCartesian3D(nx, ny, nz, Element::HEXAHEDRON,
                                      1.0, 1.0, 1.0, true);
    int dim = mesh.Dimension();

    std::cout << "Mesh: " << nx << "x" << ny << "x" << nz << " = "
              << mesh.GetNE() << " elements" << std::endl;

    // ============================================================
    // STEP 3: FINITE ELEMENT SPACES
    // ============================================================

    // Electric potential: H1 (continuous, scalar)
    H1_FECollection fec_potential(order, dim);
    FiniteElementSpace fes_potential(&mesh, &fec_potential);

    // Temperature: H1 (continuous, scalar)
    H1_FECollection fec_temp(order, dim);
    FiniteElementSpace fes_temp(&mesh, &fec_temp);

    std::cout << "Electric DOFs:   " << fes_potential.GetTrueVSize() << std::endl;
    std::cout << "Thermal DOFs:    " << fes_temp.GetTrueVSize() << std::endl;
    std::cout << "========================================\n"
              << std::endl;

    // ============================================================
    // STEP 4: BOUNDARY CONDITIONS
    // ============================================================

    // Electric BC: voltage applied at x=0 and x=L
    // Attribute 5: x=0 (ground, φ=0)
    // Attribute 3: x=L (applied voltage, φ=V)
    Array<int> elec_ess_bdr(mesh.bdr_attributes.Max());
    elec_ess_bdr = 0;
    elec_ess_bdr[4] = 1; // x=0
    elec_ess_bdr[2] = 1; // x=L

    Array<int> elec_ess_tdof_list;
    fes_potential.GetEssentialTrueDofs(elec_ess_bdr, elec_ess_tdof_list);

    // Thermal BC: convection on all boundaries (Robin BC)
    Array<int> thermal_conv_bdr(mesh.bdr_attributes.Max());
    thermal_conv_bdr = 1; // All boundaries

    // ============================================================
    // STEP 5: SOLUTION VECTORS (GridFunctions)
    // ============================================================

    GridFunction potential(&fes_potential);
    GridFunction temp_old(&fes_temp);
    GridFunction temp_new(&fes_temp);

    // Initialize
    potential = 0.0;
    temp_old = T_ref; // Start at reference temperature
    temp_new = T_ref;

    // Apply electric boundary conditions
    Array<int> bc_x0(mesh.bdr_attributes.Max());
    bc_x0 = 0;
    bc_x0[4] = 1;
    Array<int> bc_xL(mesh.bdr_attributes.Max());
    bc_xL = 0;
    bc_xL[2] = 1;

    ConstantCoefficient zero_volt(0.0);
    ConstantCoefficient applied_volt(V_applied);

    potential.ProjectBdrCoefficient(zero_volt, bc_x0);
    potential.ProjectBdrCoefficient(applied_volt, bc_xL);

    // ============================================================
    // STEP 6: OUTPUT SETUP
    // ============================================================

    ParaViewDataCollection paraview_dc("results/ex5_coupled_problem", &mesh);
    paraview_dc.SetPrefixPath("results");
    paraview_dc.SetLevelsOfDetail(order);
    paraview_dc.SetCycle(0);
    paraview_dc.SetTime(0.0);
    paraview_dc.RegisterField("Potential", &potential);
    paraview_dc.RegisterField("Temperature", &temp_new);
    paraview_dc.Save();

    // ============================================================
    // STEP 7: TIME STEPPING LOOP
    // ============================================================

    std::cout << "Starting time integration..." << std::endl;
    std::cout << "Time | Max Temp (K) | Avg Current Density (A/m²) | Coupling Iter" << std::endl;
    std::cout << "-----+-------------+---------------------------+---------------" << std::endl;

    double t = 0.0;
    int step = 0;

    while (t < t_final)
    {
        t += dt;
        step++;

        // ========================================================
        // STEP 8: COUPLING ITERATION (Picard/Fixed-point)
        // ========================================================

        temp_new = temp_old; // Initial guess

        for (int iter = 0; iter < max_coupling_iter; iter++)
        {
            // ====================================================
            // STEP 8.1: SOLVE ELECTRIC PROBLEM
            // ====================================================
            // -∇·(σ(T_new)∇φ) = 0

            // Update conductivity based on current temperature
            TemperatureDependentConductivity sigma_coef(temp_new, sigma_0, alpha_elec, T_ref);

            // Assemble electric bilinear form
            BilinearForm a_elec(&fes_potential);
            a_elec.AddDomainIntegrator(new DiffusionIntegrator(sigma_coef));
            a_elec.Assemble();
            a_elec.Finalize();

            // Form and solve linear system
            LinearForm b_elec(&fes_potential);
            b_elec = 0.0;
            b_elec.Assemble();

            SparseMatrix A_elec;
            Vector B_elec, X_elec;
            a_elec.FormLinearSystem(elec_ess_tdof_list, potential, b_elec,
                                    A_elec, X_elec, B_elec);

            GSSmoother M_elec(A_elec);
            CGSolver solver_elec;
            solver_elec.SetRelTol(1e-12);
            solver_elec.SetMaxIter(1000);
            solver_elec.SetPrintLevel(0);
            solver_elec.SetPreconditioner(M_elec);
            solver_elec.SetOperator(A_elec);
            solver_elec.Mult(B_elec, X_elec);

            a_elec.RecoverFEMSolution(X_elec, b_elec, potential);

            // ====================================================
            // STEP 8.2: COMPUTE JOULE HEATING
            // ====================================================
            // Q = σ(T)|∇φ|²

            JouleHeatingCoefficient joule_heating(potential, temp_new,
                                                  sigma_0, alpha_elec, T_ref);

            // ====================================================
            // STEP 8.3: SOLVE THERMAL PROBLEM
            // ====================================================
            // ρc(T_new - T_old)/dt + ∇·(k∇T_new) = Q + convection

            // Mass matrix: ρc/dt ∫ u v dx
            BilinearForm m_thermal(&fes_temp);
            ConstantCoefficient mass_coef(rho * cp / dt);
            m_thermal.AddDomainIntegrator(new MassIntegrator(mass_coef));
            m_thermal.Assemble();

            // Stiffness matrix: k ∫ ∇u·∇v dx
            BilinearForm stiff_thermal(&fes_temp);
            ConstantCoefficient k_coef(k_thermal);
            stiff_thermal.AddDomainIntegrator(new DiffusionIntegrator(k_coef));
            stiff_thermal.Assemble();

            // Robin BC (convection): h ∫_Γ u v ds
            ConstantCoefficient h_coef(h_conv);
            stiff_thermal.AddBoundaryIntegrator(new MassIntegrator(h_coef), thermal_conv_bdr);
            stiff_thermal.Assemble();

            // Combine: A_thermal = M + K
            SparseMatrix M_thermal, K_thermal;
            Array<int> empty_tdof; // No essential DOFs for thermal (all Robin)
            m_thermal.FormSystemMatrix(empty_tdof, M_thermal);
            stiff_thermal.FormSystemMatrix(empty_tdof, K_thermal);

            SparseMatrix *A_thermal = Add(M_thermal, K_thermal);

            // RHS: b = M*T_old + ∫ Q dx + ∫ h*T_amb ds
            LinearForm b_thermal(&fes_temp);
            b_thermal.AddDomainIntegrator(new DomainLFIntegrator(joule_heating));

            ConstantCoefficient conv_rhs(h_conv * T_amb);
            b_thermal.AddBoundaryIntegrator(new BoundaryLFIntegrator(conv_rhs),
                                            thermal_conv_bdr);
            b_thermal.Assemble();

            // Add M*T_old to RHS
            Vector rhs_thermal(fes_temp.GetTrueVSize());
            M_thermal.Mult(temp_old, rhs_thermal);
            rhs_thermal += b_thermal;

            // Solve thermal system
            Vector temp_new_vec(fes_temp.GetTrueVSize());
            temp_new_vec = temp_new; // Initial guess

            GSSmoother M_thermal_precond(*A_thermal);
            CGSolver solver_thermal;
            solver_thermal.SetRelTol(1e-12);
            solver_thermal.SetMaxIter(1000);
            solver_thermal.SetPrintLevel(0);
            solver_thermal.SetPreconditioner(M_thermal_precond);
            solver_thermal.SetOperator(*A_thermal);
            solver_thermal.Mult(rhs_thermal, temp_new_vec);

            temp_new = temp_new_vec;

            delete A_thermal;

            // ====================================================
            // STEP 8.4: CHECK COUPLING CONVERGENCE
            // ====================================================

            Vector temp_diff(temp_new_vec.Size());
            subtract(temp_new_vec, temp_old, temp_diff);
            double temp_change = temp_diff.Norml2() / temp_old.Norml2();

            if (temp_change < coupling_tol || iter == max_coupling_iter - 1)
            {
                // Output statistics
                double max_temp = temp_new.Max();

                // Compute average current density (rough estimate)
                Vector grad_phi(dim);
                double J_avg = 0.0;
                int ne = mesh.GetNE();
                for (int e = 0; e < ne; e++)
                {
                    ElementTransformation *T = mesh.GetElementTransformation(e);
                    const IntegrationRule &ir = IntRules.Get(mesh.GetElementGeometry(e), 2 * order);
                    for (int i = 0; i < ir.GetNPoints(); i++)
                    {
                        const IntegrationPoint &ip = ir.IntPoint(i);
                        T->SetIntPoint(&ip);
                        potential.GetGradient(*T, grad_phi);
                        double sigma = sigma_0 / (1.0 + alpha_elec * (temp_new.GetValue(*T, ip) - T_ref));
                        J_avg += sigma * grad_phi.Norml2() * T->Weight() * ip.weight;
                    }
                }
                J_avg /= mesh.GetElementVolume(0) * ne; // Rough average

                printf("%4d | %11.2f | %25.2e | %13d\n",
                       step, max_temp, J_avg, iter + 1);
                break;
            }
        }

        // Update for next time step
        temp_old = temp_new;

        // Save output every 10 steps
        if (step % 10 == 0)
        {
            paraview_dc.SetCycle(step);
            paraview_dc.SetTime(t);
            paraview_dc.Save();
        }
    }

    std::cout << "========================================" << std::endl;
    std::cout << "Simulation complete!" << std::endl;
    std::cout << "Output: results/ex5_coupled_problem/" << std::endl;
    std::cout << "========================================" << std::endl;

    return 0;
}
