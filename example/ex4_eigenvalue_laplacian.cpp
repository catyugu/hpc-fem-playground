#include "mfem.hpp"
#include <fstream>
#include <iostream>

// Eigenvalue problem example: Laplacian eigenvalue problem
// Problem: -∇²u = λu in Ω
//          u = 0 on ∂Ω (homogeneous Dirichlet BC)
//
// Weak form: ∫ ∇u·∇v dx = λ ∫ u v dx
// Discrete form: K u = λ M u
//
// This computes the smallest eigenvalues (vibration modes, quantum states, etc.)
// using LOBPCG (Locally Optimal Block Preconditioned Conjugate Gradient)

using namespace mfem;

int main(int argc, char *argv[])
{
    // --- (I) Initialize MPI if available ------------------------------------
    int num_procs = 1;
    int myid = 0;

#ifdef MFEM_USE_MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
#endif

    // --- (II) Command line options ------------------------------------------
    const char *mesh_file = "testdata/testmesh_cube.mesh";
    const char *output_dir = "results/ex4_eigenvalue";
    int order = 1;
    int nev = 5;   // Number of eigenvalues to compute
    int seed = 75; // Random seed for LOBPCG initial vectors

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
    if (argc > 4)
    {
        nev = atoi(argv[4]);
    }
    if (argc > 5)
    {
        seed = atoi(argv[5]);
    }

    if (myid == 0)
    {
        std::cout << "========================================" << std::endl;
        std::cout << "Laplacian Eigenvalue Problem Example" << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Mesh file:       " << mesh_file << std::endl;
        std::cout << "Output dir:      " << output_dir << std::endl;
        std::cout << "FE order:        " << order << std::endl;
        std::cout << "Num eigenvalues: " << nev << std::endl;
        std::cout << "Random seed:     " << seed << std::endl;
        std::cout << "========================================" << std::endl;
    }

    // --- (III) Mesh and FE space --------------------------------------------
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    if (myid == 0)
    {
        std::cout << "Mesh dimension:  " << dim << std::endl;
        std::cout << "Mesh elements:   " << mesh->GetNE() << std::endl;
    }

    // Define finite element space (H1 for Laplacian)
    H1_FECollection fec(order, dim);
    FiniteElementSpace *fes = new FiniteElementSpace(mesh, &fec);

    if (myid == 0)
    {
        std::cout << "Number of DOFs:  " << fes->GetTrueVSize() << std::endl;
    }

    // --- (IV) Bilinear forms: Stiffness and Mass ----------------------------
    // Stiffness matrix K: ∫ ∇u·∇v dx
    BilinearForm *K = new BilinearForm(fes);
    ConstantCoefficient one(1.0);
    K->AddDomainIntegrator(new DiffusionIntegrator(one));
    K->Assemble();
    K->Finalize();

    // Mass matrix M: ∫ u v dx
    BilinearForm *M = new BilinearForm(fes);
    M->AddDomainIntegrator(new MassIntegrator(one));
    M->Assemble();
    M->Finalize();

    // --- (V) Essential (Dirichlet) boundary conditions ----------------------
    // Homogeneous Dirichlet on all boundaries: u = 0
    Array<int> ess_bdr(mesh->bdr_attributes.Max());
    ess_bdr = 1; // Mark all boundaries as essential

    Array<int> ess_tdof_list;
    fes->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);

    if (myid == 0)
    {
        std::cout << "Essential DOFs:  " << ess_tdof_list.Size() << std::endl;
    }

    // Form constrained matrices (eliminate boundary DOFs)
    SparseMatrix K_mat, M_mat;
    K->FormSystemMatrix(ess_tdof_list, K_mat);
    M->FormSystemMatrix(ess_tdof_list, M_mat);

    // --- (VI) Setup eigensolver ---------------------------------------------
#ifdef MFEM_USE_MPI
    // Parallel version using HypreLOBPCG
    if (myid == 0)
    {
        std::cout << "\nSetting up LOBPCG eigensolver (MPI version)..." << std::endl;
    }

    // Convert to HypreParMatrix (parallel format)
    // For serial run as MPI with 1 process
    HYPRE_BigInt glob_size = fes->GetTrueVSize();
    HYPRE_BigInt row_starts[2] = {0, glob_size};

    HypreParMatrix *K_hypre = new HypreParMatrix(MPI_COMM_WORLD,
                                                 glob_size, row_starts,
                                                 &K_mat);
    HypreParMatrix *M_hypre = new HypreParMatrix(MPI_COMM_WORLD,
                                                 glob_size, row_starts,
                                                 &M_mat);

    // Setup preconditioner (BoomerAMG)
    HypreBoomerAMG *amg = new HypreBoomerAMG(*K_hypre);
    amg->SetPrintLevel(0);

    // LOBPCG eigensolver
    HypreLOBPCG *lobpcg = new HypreLOBPCG(MPI_COMM_WORLD);
    lobpcg->SetNumModes(nev);
    lobpcg->SetPreconditioner(*amg);
    lobpcg->SetMaxIter(200);
    lobpcg->SetTol(1e-8);
    lobpcg->SetPrecondUsageMode(1);
    lobpcg->SetPrintLevel(1);
    lobpcg->SetMassMatrix(*M_hypre);
    lobpcg->SetOperator(*K_hypre);

    // Set initial random vectors
    srand(seed);
    HypreParVector **init_vecs = new HypreParVector *[nev];
    for (int i = 0; i < nev; i++)
    {
        init_vecs[i] = new HypreParVector(MPI_COMM_WORLD, glob_size, row_starts);
        init_vecs[i]->Randomize(seed + i);
    }
    lobpcg->SetInitialVectors(nev, init_vecs);

    // --- (VII) Solve eigenvalue problem -------------------------------------
    if (myid == 0)
    {
        std::cout << "\nSolving eigenvalue problem..." << std::endl;
    }

    lobpcg->Solve();

    // --- (VIII) Extract eigenvalues -----------------------------------------
    Array<double> eigenvalues;
    lobpcg->GetEigenvalues(eigenvalues);

    if (myid == 0)
    {
        std::cout << "\nComputed eigenvalues:" << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        std::cout << "  Mode  |  Eigenvalue  |  sqrt(λ)  |  f=√λ/(2π)" << std::endl;
        std::cout << "---------------------------------------------" << std::endl;
        for (int i = 0; i < nev && i < eigenvalues.Size(); i++)
        {
            double lambda = eigenvalues[i];
            double sqrt_lambda = sqrt(lambda);
            double freq = sqrt_lambda / (2.0 * M_PI);
            std::cout << std::setw(6) << i << "  | "
                      << std::setw(12) << std::scientific << std::setprecision(5) << lambda << " | "
                      << std::setw(9) << std::fixed << std::setprecision(4) << sqrt_lambda << " | "
                      << std::setw(9) << std::setprecision(4) << freq << std::endl;
        }
        std::cout << "---------------------------------------------" << std::endl;
    }

    // --- (IX) Extract and save eigenvectors ---------------------------------
    if (myid == 0)
    {
        std::cout << "\nSaving eigenmodes to ParaView..." << std::endl;
    }

    ParaViewDataCollection paraview(output_dir, mesh);
    paraview.SetLevelsOfDetail(order);

    for (int i = 0; i < nev && i < eigenvalues.Size(); i++)
    {
        GridFunction u(fes);
        u = lobpcg->GetEigenvector(i);

        // Normalize
        double norm = u.Norml2();
        if (norm > 0.0)
        {
            u /= norm;
        }

        // Create field name
        std::ostringstream mode_name;
        mode_name << "mode_" << std::setw(3) << std::setfill('0') << i;

        // Register and save (use eigenvalue as "time" for easy identification)
        paraview.RegisterField(mode_name.str().c_str(), &u);
        paraview.SetCycle(i);
        paraview.SetTime(eigenvalues[i]);
        paraview.Save();
    }

    // Cleanup
    for (int i = 0; i < nev; i++)
    {
        delete init_vecs[i];
    }
    delete[] init_vecs;
    delete lobpcg;
    delete amg;
    delete K_hypre;
    delete M_hypre;

#else
    // Sequential version (no MPI)
    if (myid == 0)
    {
        std::cout << "\nMFEM not compiled with MPI support." << std::endl;
        std::cout << "Cannot use HypreLOBPCG for eigenvalue problems." << std::endl;
        std::cout << "Please recompile MFEM with HYPRE and MPI enabled." << std::endl;
    }
#endif

    // --- (X) Validation against analytical solution (for unit cube) ---------
#ifdef MFEM_USE_MPI
    if (myid == 0 && dim == 3)
    {
        // For unit cube [0,1]^3 with Dirichlet BC, analytical eigenvalues are:
        // λ_{l,m,n} = π²(l² + m² + n²) for l,m,n = 1,2,3,...

        std::cout << "\n=== Validation (Unit Cube [0,1]^3) ===" << std::endl;
        std::cout << "Analytical eigenvalues: λ = π²(l² + m² + n²)" << std::endl;
        std::cout << "\nExpected first few modes:" << std::endl;
        std::cout << "  (1,1,1): λ ≈ " << M_PI * M_PI * 3.0 << std::endl;
        std::cout << "  (1,1,2): λ ≈ " << M_PI * M_PI * 6.0 << std::endl;
        std::cout << "  (1,2,2): λ ≈ " << M_PI * M_PI * 9.0 << std::endl;

        // Check mesh bounds
        Vector bb_min, bb_max;
        mesh->GetBoundingBox(bb_min, bb_max);
        bool is_unit_cube = true;
        double tol = 1e-3;
        for (int i = 0; i < dim; i++)
        {
            if (fabs(bb_min(i) - 0.0) > tol || fabs(bb_max(i) - 1.0) > tol)
            {
                is_unit_cube = false;
                break;
            }
        }

        if (is_unit_cube)
        {
            std::cout << "\nMesh appears to be unit cube." << std::endl;
            std::cout << "Relative errors:" << std::endl;

            std::vector<double> analytical = {
                M_PI * M_PI * 3.0, // (1,1,1)
                M_PI * M_PI * 6.0, // (1,1,2), (1,2,1), (2,1,1) - degenerate
                M_PI * M_PI * 9.0, // (1,2,2), (2,1,2), (2,2,1) - degenerate
            };

            for (size_t i = 0; i < analytical.size() && i < (size_t)eigenvalues.Size(); i++)
            {
                double err = fabs(eigenvalues[i] - analytical[i]) / analytical[i];
                std::cout << "  Mode " << i << ": "
                          << std::setw(10) << std::scientific << err;
                if (err < 0.01)
                {
                    std::cout << " (good)";
                }
                else if (err < 0.1)
                {
                    std::cout << " (acceptable)";
                }
                else
                {
                    std::cout << " (check mesh/refinement)";
                }
                std::cout << std::endl;
            }
        }
        else
        {
            std::cout << "\nMesh is not unit cube. Skipping analytical validation." << std::endl;
        }
    }
#endif

    // --- (XI) Finalize ------------------------------------------------------
    if (myid == 0)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "Eigenvalue computation completed!" << std::endl;
        std::cout << "Output saved to: " << output_dir << std::endl;
        std::cout << "========================================" << std::endl;
    }

    delete M;
    delete K;
    delete fes;
    delete mesh;

#ifdef MFEM_USE_MPI
    MPI_Finalize();
#endif

    return 0;
}
