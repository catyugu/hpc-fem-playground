#include "mfem.hpp"
#include <fstream>
#include <iostream>

// Transient (time-dependent) heat equation example
// Problem: du/dt - div(k grad u) = 0 in Omega x (0, T]
//          u(x, 0) = u0(x)
//          u = g on Gamma_D (Dirichlet BC)
//          -k grad u · n = q on Gamma_N (Neumann BC)
//
// Weak form: M du/dt + K u = b
// where M = mass matrix, K = stiffness matrix, b = RHS (boundary fluxes)
//
// Time integration: Backward Euler (implicit) or RK4 (explicit)

using namespace mfem;

// Initial condition function
double initial_condition(const Vector &x)
{
    // Gaussian temperature distribution
    double x0 = 0.5, y0 = 0.5, z0 = 0.5;  // Center
    double sigma = 0.1;  // Width
    
    double r2 = 0.0;
    r2 += (x.Size() >= 1) ? (x(0) - x0) * (x(0) - x0) : 0.0;
    r2 += (x.Size() >= 2) ? (x(1) - y0) * (x(1) - y0) : 0.0;
    r2 += (x.Size() >= 3) ? (x(2) - z0) * (x(2) - z0) : 0.0;
    
    return 500.0 * exp(-r2 / (2.0 * sigma * sigma)) + 293.15;
}

// Time-dependent operator for heat equation: M du/dt = -K u + b
class HeatOperator : public TimeDependentOperator
{
private:
    BilinearForm &M_form;      // Mass matrix form
    BilinearForm &K_form;      // Stiffness matrix form
    const Vector &b;           // RHS vector (boundary terms)
    
    CGSolver M_solver;         // Mass matrix solver (for explicit)
    mutable Vector z;          // Auxiliary vector
    mutable SparseMatrix *A_mat;  // System matrix for implicit (M + dt*K)
    mutable double dt_cached;     // Cached time step
    
public:
    HeatOperator(BilinearForm &M_, BilinearForm &K_, const Vector &b_)
        : TimeDependentOperator(M_.Size()), M_form(M_), K_form(K_), b(b_),
          A_mat(nullptr), dt_cached(-1.0)
    {
        z.SetSize(M_form.Size());
        
        // Setup mass matrix solver for explicit methods
        M_solver.SetOperator(M_form.SpMat());
        M_solver.SetRelTol(1e-9);
        M_solver.SetMaxIter(100);
        M_solver.SetPrintLevel(0);
    }
    
    // Explicit time stepping: du/dt = M^{-1}(-K*u + b)
    virtual void Mult(const Vector &u, Vector &dudt) const override
    {
        K_form.SpMat().Mult(u, z);  // z = K*u
        z.Neg();                     // z = -K*u
        z += b;                      // z = -K*u + b
        M_solver.Mult(z, dudt);      // dudt = M^{-1} * z
    }
    
    // Implicit solve: given dt, u_old, solve for k where u_new = u_old + dt*k
    // Need to solve: (M + dt*K)*k = -K*u_old + b
    virtual void ImplicitSolve(const double dt, const Vector &u, Vector &k) override
    {
        // Form system matrix if dt changed: A = M + dt*K
        if (dt != dt_cached) {
            delete A_mat;
            A_mat = Add(1.0, M_form.SpMat(), dt, K_form.SpMat());
            dt_cached = dt;
        }
        
        // RHS: rhs = -K*u + b
        Vector rhs(u.Size());
        K_form.SpMat().Mult(u, rhs);
        rhs.Neg();
        rhs += b;
        
        // Solve A*k = rhs
        CGSolver solver;
        solver.SetOperator(*A_mat);
        solver.SetRelTol(1e-9);
        solver.SetMaxIter(500);
        solver.SetPrintLevel(0);
        solver.Mult(rhs, k);
    }
    
    virtual ~HeatOperator() { delete A_mat; }
};

int main(int argc, char *argv[])
{
    // --- (I) Command line options -------------------------------------------
    const char *mesh_file = "testdata/testmesh_cube.mesh";
    const char *output_dir = "results/ex3_transient_heat";
    int order = 1;
    double t_final = 1.0;
    double dt = 0.01;
    int ode_solver_type = 2;  // 1 = Forward Euler, 2 = Backward Euler, 3 = RK4
    int vis_steps = 10;        // Output every N steps
    
    if (argc > 1) { mesh_file = argv[1]; }
    if (argc > 2) { output_dir = argv[2]; }
    if (argc > 3) { order = atoi(argv[3]); }
    if (argc > 4) { t_final = atof(argv[4]); }
    if (argc > 5) { dt = atof(argv[5]); }
    if (argc > 6) { ode_solver_type = atoi(argv[6]); }
    
    std::cout << "========================================" << std::endl;
    std::cout << "Transient Heat Equation Example" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Mesh file:     " << mesh_file << std::endl;
    std::cout << "Output dir:    " << output_dir << std::endl;
    std::cout << "FE order:      " << order << std::endl;
    std::cout << "Final time:    " << t_final << std::endl;
    std::cout << "Time step:     " << dt << std::endl;
    std::cout << "ODE solver:    " << ode_solver_type 
              << " (1=FwdEuler, 2=BwdEuler, 3=RK4)" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // --- (II) Mesh and FE space ---------------------------------------------
    Mesh mesh(mesh_file);
    int dim = mesh.Dimension();
    std::cout << "Mesh dimension: " << dim << std::endl;
    std::cout << "Mesh elements:  " << mesh.GetNE() << std::endl;
    
    H1_FECollection fec(order, dim);
    FiniteElementSpace fes(&mesh, &fec);
    std::cout << "Number of DOFs: " << fes.GetTrueVSize() << std::endl;
    
    // --- (III) Material properties ------------------------------------------
    double k = 45.0;  // Thermal conductivity [W/m/K] (e.g., steel)
    double rho = 7850.0;  // Density [kg/m^3]
    double cp = 460.0;    // Specific heat [J/kg/K]
    double alpha = k / (rho * cp);  // Thermal diffusivity [m^2/s]
    
    std::cout << "Thermal diffusivity: " << alpha << " m^2/s" << std::endl;
    
    // --- (IV) Bilinear forms: Mass and Stiffness ----------------------------
    BilinearForm M(&fes);
    ConstantCoefficient one(1.0);
    M.AddDomainIntegrator(new MassIntegrator(one));
    M.Assemble();
    M.Finalize();
    
    BilinearForm K(&fes);
    ConstantCoefficient k_coef(k);
    K.AddDomainIntegrator(new DiffusionIntegrator(k_coef));
    K.Assemble();
    K.Finalize();
    
    // --- (V) Linear form (RHS): boundary fluxes -----------------------------
    LinearForm b(&fes);
    b = 0.0;  // Initialize to zero
    // No boundary flux for this simple example
    b.Assemble();
    
    // --- (VI) Initial condition ---------------------------------------------
    FunctionCoefficient u0_func(initial_condition);
    GridFunction u(&fes);
    u.ProjectCoefficient(u0_func);
    
    std::cout << "Initial temperature: min = " << u.Min() 
              << ", max = " << u.Max() << std::endl;
    
    // --- (VII) Essential (Dirichlet) BCs ------------------------------------
    // Example: Fix temperature on attribute 1 to ambient
    Array<int> ess_bdr(mesh.bdr_attributes.Max());
    ess_bdr = 0;
    
    std::cout << "Boundary attributes: " << mesh.bdr_attributes.Max() << std::endl;
    for (int i = 0; i < mesh.bdr_attributes.Size(); i++) {
        std::cout << "  Attr " << mesh.bdr_attributes[i] << std::endl;
    }
    
    if (mesh.bdr_attributes.Max() >= 1) {
        ess_bdr[1-1] = 1;  // Attribute 1: fixed temperature
    }
    
    Array<int> ess_tdof_list;
    fes.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
    
    std::cout << "Essential DOFs: " << ess_tdof_list.Size() << std::endl;
    
    // Project Dirichlet values
    ConstantCoefficient ambient_temp(293.15);  // 20°C
    if (ess_tdof_list.Size() > 0) {
        u.ProjectBdrCoefficient(ambient_temp, ess_bdr);
        std::cout << "Applied Dirichlet BC on " << ess_tdof_list.Size() 
                  << " DOFs" << std::endl;
    }
    
    // --- (VIII) Time-dependent operator -------------------------------------
    HeatOperator oper(M, K, b);
    
    // --- (IX) ODE solver setup ----------------------------------------------
    ODESolver *ode_solver = nullptr;
    switch (ode_solver_type) {
        case 1:
            ode_solver = new ForwardEulerSolver;
            std::cout << "Using Forward Euler (explicit)" << std::endl;
            break;
        case 2:
            ode_solver = new BackwardEulerSolver;
            std::cout << "Using Backward Euler (implicit)" << std::endl;
            break;
        case 3:
            ode_solver = new RK4Solver;
            std::cout << "Using RK4 (explicit)" << std::endl;
            break;
        default:
            std::cerr << "Unknown ODE solver type: " << ode_solver_type << std::endl;
            return 1;
    }
    
    ode_solver->Init(oper);
    
    // --- (X) Time integration loop ------------------------------------------
    double t = 0.0;
    bool done = false;
    int ti = 0;
    
    // Setup ParaView output
    ParaViewDataCollection paraview(output_dir, &mesh);
    paraview.SetLevelsOfDetail(order);
    paraview.RegisterField("Temperature", &u);
    paraview.SetTime(t);
    paraview.SetCycle(ti);
    paraview.Save();
    
    std::cout << "\nTime integration started..." << std::endl;
    std::cout << "Step |   Time    | T_min  | T_max  | T_avg" << std::endl;
    std::cout << "-----+-----------+--------+--------+--------" << std::endl;
    
    while (!done) {
        // Adjust time step for final step
        if (t + dt >= t_final - dt/2) {
            dt = t_final - t;
            done = true;
        }
        
        // Take time step
        ode_solver->Step(u, t, dt);
        ti++;
        
        // Re-apply Dirichlet BCs (important for time-dependent or implicit methods)
        if (ess_tdof_list.Size() > 0) {
            u.ProjectBdrCoefficient(ambient_temp, ess_bdr);
        }
        
        // Output
        if (ti % vis_steps == 0 || done) {
            double u_min = u.Min();
            double u_max = u.Max();
            double u_avg = u.Sum() / u.Size();
            
            std::cout << std::setw(4) << ti << " | "
                      << std::setw(9) << std::fixed << std::setprecision(4) << t << " | "
                      << std::setw(6) << std::setprecision(1) << u_min << " | "
                      << std::setw(6) << std::setprecision(1) << u_max << " | "
                      << std::setw(6) << std::setprecision(1) << u_avg << std::endl;
            
            paraview.SetTime(t);
            paraview.SetCycle(ti);
            paraview.Save();
        }
    }
    
    std::cout << "\nTime integration completed!" << std::endl;
    std::cout << "Total steps: " << ti << std::endl;
    std::cout << "Final time:  " << t << std::endl;
    std::cout << "Output saved to: " << output_dir << std::endl;
    
    // --- (XI) Cleanup -------------------------------------------------------
    delete ode_solver;
    
    return 0;
}
