/**
 * @file ex8_maxwell_rk4.cpp
 * @brief Example 8: Time-domain Maxwell solver with RK4 integration
 * 
 * This example demonstrates the time-domain Maxwell solver using MFEM's
 * 4th-order Runge-Kutta time integrator for improved energy conservation.
 * 
 * The coupled Maxwell equations are:
 *   ε ∂E/∂t = ∇×(μ⁻¹B) - σE - J
 *   ∂B/∂t = -∇×E
 * 
 * We solve this as a system dy/dt = f(t,y) where y = [E; B].
 * 
 * Compile: make ex8_maxwell_rk4
 * Run (serial): ./ex8_maxwell_rk4
 * Run (parallel): mpirun -np 4 ./ex8_maxwell_rk4
 */

#include "mfem.hpp"
#include "hpcfem/physics/physics_maxwell_timedomain.hpp"
#include <iostream>
#include <fstream>

using namespace std;
using namespace mfem;
using namespace hpcfem;

// Material properties (free space)
double epsilon(const Vector& x)
{
    return 8.8541878176e-12;  // F/m
}

double muInv(const Vector& x)
{
    return 1.0 / (4.0e-7 * M_PI);  // H^-1 m^-1
}

// Initial E-field: simple cavity mode
class InitialEField : public VectorCoefficient
{
public:
    InitialEField() : VectorCoefficient(3) {}
    
    virtual void Eval(Vector& V, ElementTransformation& T,
                     const IntegrationPoint& ip)
    {
        double x[3];
        Vector transip(x, 3);
        T.Transform(ip, transip);
        
        V.SetSize(3);
        V(0) = 0.0;
        V(1) = 0.0;
        // E_z = sin(pi*x) * sin(pi*y)
        V(2) = sin(M_PI * x[0]) * sin(M_PI * x[1]);
    }
};

/**
 * @brief Wrapper operator for coupled Maxwell system
 * 
 * Combines both E and B fields into a single state vector for
 * MFEM's ODE solvers. The state vector is organized as [E; B].
 */
class MaxwellODEOperator : public TimeDependentOperator
{
public:
    MaxwellODEOperator(PhysicsMaxwellTimeDomain* solver)
        : TimeDependentOperator(solver->getEVector()->Size() + 
                               solver->getBVector()->Size()),
          solver_(solver),
          eSize_(solver->getEVector()->Size()),
          bSize_(solver->getBVector()->Size()),
          dEdt_(eSize_),
          dBdt_(bSize_)
    {
        // Total state size = E DOFs + B DOFs
    }
    
    virtual void Mult(const Vector& y, Vector& dydt) const override
    {
        // Extract E and B from state vector
        Vector E(y.GetData(), eSize_);
        Vector B(y.GetData() + eSize_, bSize_);
        
        // Get output views
        Vector dEdt_view(dydt.GetData(), eSize_);
        Vector dBdt_view(dydt.GetData() + eSize_, bSize_);
        
        // Compute dE/dt = M_ε^{-1} [∇×(μ⁻¹B) - σE - J]
        solver_->Mult(B, dEdt_);
        dEdt_view = dEdt_;
        
        // Compute dB/dt = -∇×E
        solver_->getNegCurlOperator()->Mult(E, dBdt_);
        dBdt_view = dBdt_;
    }
    
    void GetState(Vector& state) const
    {
        // Copy current E and B into state vector
        Vector* E = solver_->getEVector();
        Vector* B = solver_->getBVector();
        
        Vector E_view(state.GetData(), eSize_);
        Vector B_view(state.GetData() + eSize_, bSize_);
        
        E_view = *E;
        B_view = *B;
    }
    
    void SetState(const Vector& state)
    {
        // Copy state vector into E and B
        Vector* E = solver_->getEVector();
        Vector* B = solver_->getBVector();
        
        const Vector E_view(state.GetData(), eSize_);
        const Vector B_view(state.GetData() + eSize_, bSize_);
        
        *E = E_view;
        *B = B_view;
        
        solver_->syncGridFunctions();
    }
    
private:
    PhysicsMaxwellTimeDomain* solver_;
    int eSize_;
    int bSize_;
    mutable Vector dEdt_;
    mutable Vector dBdt_;
};

int main(int argc, char *argv[])
{
#ifdef MFEM_USE_MPI
    MPI_Init(&argc, &argv);
    Hypre::Init();
#endif

    // Parse command line options
    int order = 1;
    int refLevels = 2;
    double tFinal = 1.0e-9;  // 1 nanosecond
    int numSteps = 20;
    bool visualization = false;
    
    OptionsParser args(argc, argv);
    args.AddOption(&order, "-o", "--order",
                   "Finite element order (polynomial degree).");
    args.AddOption(&refLevels, "-r", "--refine",
                   "Number of times to refine the mesh uniformly.");
    args.AddOption(&tFinal, "-tf", "--t-final",
                   "Final time; start time is 0.");
    args.AddOption(&numSteps, "-n", "--num-steps",
                   "Number of time steps.");
    args.AddOption(&visualization, "-vis", "--visualization", "-no-vis",
                   "--no-visualization",
                   "Enable or disable GLVis visualization.");
    args.Parse();
    if (!args.Good())
    {
        args.PrintUsage(cout);
        return 1;
    }
    args.PrintOptions(cout);

    int myRank = 0;
    int numProcs = 1;
    
#ifdef MFEM_USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
#endif

    // Create mesh: unit cube cavity
    Mesh* mesh = new Mesh(Mesh::MakeCartesian3D(
        4, 4, 4, Element::HEXAHEDRON, 1.0, 1.0, 1.0));
    
    for (int i = 0; i < refLevels; i++)
    {
        mesh->UniformRefinement();
    }

#ifdef MFEM_USE_MPI
    ParMesh* pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
#endif

    if (myRank == 0)
    {
        cout << "\nCreating Maxwell solver with RK4 integrator..." << endl;
    }

    // Create Maxwell solver (lossless cavity)
    Array<int> emptyArray;  // No special boundary conditions
    
#ifdef MFEM_USE_MPI
    PhysicsMaxwellTimeDomain solver(pmesh, order,
                                     epsilon, muInv,
                                     nullptr,  // No losses
                                     nullptr,  // No sources
                                     emptyArray, emptyArray,
                                     nullptr);
#else
    PhysicsMaxwellTimeDomain solver(mesh, order,
                                     epsilon, muInv,
                                     nullptr, nullptr,
                                     emptyArray, emptyArray,
                                     nullptr);
#endif

    // Set initial conditions
    InitialEField e0;
    solver.setInitialEField(e0);

    // Create ODE operator wrapper
    MaxwellODEOperator odeOp(&solver);
    
    // Get initial state
    Vector state(odeOp.Height());
    odeOp.GetState(state);
    
    double initialEnergy = solver.getEnergy();
    
    if (myRank == 0)
    {
        cout << "\nInitial energy: " << initialEnergy << " J" << endl;
        cout << "Using 4th-order Runge-Kutta time integration" << endl;
    }

    // Create RK4 time stepper
    RK4Solver rkSolver;
    rkSolver.Init(odeOp);
    
    // Time step size
    double dt = tFinal / numSteps;
    double t = 0.0;
    
    if (myRank == 0)
    {
        cout << "Time step: " << dt << " s" << endl;
        cout << "Number of steps: " << numSteps << endl;
        cout << "\nTime stepping with RK4..." << endl;
    }
    
    double maxEnergyError = 0.0;
    
    // Time stepping loop
    for (int step = 0; step < numSteps; step++)
    {
        // Take one RK4 step
        rkSolver.Step(state, t, dt);
        
        // Update solver state
        odeOp.SetState(state);
        
        // Update source time if needed
        solver.updateSourceTime(t);
        
        // Compute energy
        double energy = solver.getEnergy();
        double energyError = fabs(energy - initialEnergy) / initialEnergy;
        
        if (energyError > maxEnergyError)
        {
            maxEnergyError = energyError;
        }
        
        if (myRank == 0 && numSteps > 0 && step % max(1, numSteps / 10) == 0)
        {
            cout << "Step " << step << "/" << numSteps 
                 << ", t = " << t << " s"
                 << ", Energy = " << energy << " J"
                 << ", Error = " << energyError << endl;
        }
    }

    // Get final energy (collective operation)
    double finalEnergy = solver.getEnergy();

    if (myRank == 0)
    {
        cout << "\nSimulation complete!" << endl;
        cout << "Final energy: " << finalEnergy << " J" << endl;
        cout << "Initial energy: " << initialEnergy << " J" << endl;
        cout << "Max energy error: " << maxEnergyError << endl;
        
        if (maxEnergyError < 1e-6)
        {
            cout << "✓ Energy conserved to high precision!" << endl;
        }
        else if (maxEnergyError < 1e-3)
        {
            cout << "✓ Energy well conserved (< 0.1%)" << endl;
        }
        else if (maxEnergyError < 0.01)
        {
            cout << "✓ Energy reasonably conserved (< 1%)" << endl;
        }
        else
        {
            cout << "⚠ Energy conservation could be improved" << endl;
        }
    }

#ifdef MFEM_USE_MPI
    delete pmesh;
    MPI_Finalize();
#else
    delete mesh;
#endif

    return 0;
}
