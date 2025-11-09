/**
 * @file ex7_maxwell_cavity.cpp
 * @brief Example 7: Time-domain Maxwell solver - Cavity resonance
 * 
 * This example demonstrates the time-domain Maxwell solver with energy
 * conservation in a lossless PEC cavity.
 * 
 * Compile: make ex7_maxwell_cavity
 * Run (serial): ./ex7_maxwell_cavity
 * Run (parallel): mpirun -np 4 ./ex7_maxwell_cavity
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

// Initial E-field: simple mode
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
    int numSteps = 100;
    bool visualization = true;
    
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

    // Create mesh
    Mesh* mesh = new Mesh(Mesh::MakeCartesian3D(
        4, 4, 4, Element::HEXAHEDRON, 1.0, 1.0, 1.0));
    
    for (int i = 0; i < refLevels; i++)
    {
        mesh->UniformRefinement();
    }

#ifdef MFEM_USE_MPI
    ParMesh* pmesh = new ParMesh(MPI_COMM_WORLD, *mesh);
    delete mesh;
#else
    ParMesh* pmesh = mesh;
#endif

    if (myRank == 0)
    {
        cout << "Creating Maxwell solver..." << endl;
    }

    // Create solver
    Array<int> emptyArray;
    PhysicsMaxwellTimeDomain solver(pmesh, order,
                                     epsilon,
                                     muInv,
                                     nullptr,  // No conductivity (lossless)
                                     nullptr,  // No source
                                     emptyArray, emptyArray,
                                     nullptr);

    // Set initial E-field
    InitialEField e0;
    solver.setInitialEField(e0);

    // Get maximum stable time step
    double dtMax = solver.getMaximumTimeStep();
    double dt = 0.5 * dtMax;  // Safety factor
    int actualSteps = static_cast<int>(tFinal / dt);
    
    if (actualSteps < numSteps)
    {
        // User requested more steps than needed for stability
        numSteps = actualSteps;
    }
    else
    {
        // Adjust dt to match user's requested number of steps
        dt = tFinal / numSteps;
    }

    if (myRank == 0)
    {
        cout << "Max stable dt: " << dtMax << " s" << endl;
        cout << "Using dt: " << dt << " s" << endl;
        cout << "Number of steps: " << numSteps << endl;
    }

    // Get field vectors
    Vector* E = solver.getEVector();
    Vector* B = solver.getBVector();
    Vector dBdt(B->Size());
    Vector dEdt(E->Size());

    // Initial energy
    double initialEnergy = solver.getEnergy();
    
    if (myRank == 0)
    {
        cout << "\nInitial energy: " << initialEnergy << " J" << endl;
        cout << "\nTime stepping..." << endl;
    }

    // Time stepping loop
    double t = 0.0;
    double maxEnergyError = 0.0;
    
    for (int step = 0; step < numSteps; step++)
    {
        // Update B: B^{n+1} = B^n + dt * (-curl(E^n))
        solver.getNegCurlOperator()->Mult(*E, dBdt);
        B->Add(dt, dBdt);
        
        // Update E: E^{n+1} = E^n + dt * dE/dt
        solver.Mult(*B, dEdt);
        E->Add(dt, dEdt);
        
        t += dt;
        
        // Sync grid functions for energy computation
        solver.syncGridFunctions();
        
        // Check energy conservation
        double energy = solver.getEnergy();
        double energyError = abs((energy - initialEnergy) / initialEnergy);
        
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

    // Get final energy (collective operation, all ranks participate)
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
