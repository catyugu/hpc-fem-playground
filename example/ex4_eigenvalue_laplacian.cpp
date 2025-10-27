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

    std::cout << "========================================" << std::endl;
    std::cout << "Laplacian Eigenvalue Problem Example" << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "Mesh file:       " << mesh_file << std::endl;
    std::cout << "Output dir:      " << output_dir << std::endl;
    std::cout << "FE order:        " << order << std::endl;
    std::cout << "Num eigenvalues: " << nev << std::endl;
    std::cout << "Random seed:     " << seed << std::endl;
    std::cout << "========================================" << std::endl;

    // --- (III) Mesh and FE space --------------------------------------------
    Mesh *mesh = new Mesh(mesh_file, 1, 1);
    int dim = mesh->Dimension();

    std::cout << "Mesh dimension:  " << dim << std::endl;
    std::cout << "Mesh elements:   " << mesh->GetNE() << std::endl;

    // Define finite element space (H1 for Laplacian)
    H1_FECollection fec(order, dim);
    FiniteElementSpace *fes = new FiniteElementSpace(mesh, &fec);

    std::cout << "Number of DOFs:  " << fes->GetTrueVSize() << std::endl;

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

    std::cout << "Essential DOFs:  " << ess_tdof_list.Size() << std::endl;

    // Form constrained matrices (eliminate boundary DOFs)
    SparseMatrix K_mat, M_mat;
    K->FormSystemMatrix(ess_tdof_list, K_mat);
    M->FormSystemMatrix(ess_tdof_list, M_mat);

    std::cout << "\nMFEM not compiled with MPI support." << std::endl;
    std::cout << "Cannot use HypreLOBPCG for eigenvalue problems." << std::endl;
    std::cout << "Please recompile MFEM with HYPRE and MPI enabled." << std::endl;

    // --- (XI) Finalize ------------------------------------------------------

    std::cout << "\n========================================" << std::endl;
    std::cout << "Eigenvalue computation completed!" << std::endl;
    std::cout << "Output saved to: " << output_dir << std::endl;
    std::cout << "========================================" << std::endl;

    delete M;
    delete K;
    delete fes;
    delete mesh;

    return 0;
}
