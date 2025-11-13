/**
 * @file physics/physics_waveguide_eigen.cpp
 * @brief Implementation of 2D waveguide eigenmode solver (TE/TM)
 */

#include "physics_waveguide_eigen.hpp"
#include <stdexcept>
#include <iostream>

namespace hpcfem
{

#ifdef MFEM_USE_MPI
PhysicsWaveguideEigen::PhysicsWaveguideEigen(mfem::ParMesh *pmesh, int order, WaveguideModeType mode)
    : pmesh_(pmesh), stiffForm_(nullptr), massForm_(nullptr), fespace_(nullptr), fec_(nullptr), modeType_(mode), order_(order)
{
    if (!pmesh_) { throw std::runtime_error("PhysicsWaveguideEigen: null ParMesh"); }
    int dim = pmesh_->Dimension();
    if (dim != 2) { throw std::runtime_error("PhysicsWaveguideEigen: only 2D meshes supported"); }

    fec_ = new mfem::H1_FECollection(order_, dim);
    fespace_ = new mfem::ParFiniteElementSpace(pmesh_, fec_);
}
#else
PhysicsWaveguideEigen::PhysicsWaveguideEigen(mfem::Mesh *mesh, int order, WaveguideModeType mode)
    : mesh_(mesh), stiffForm_(nullptr), massForm_(nullptr), fespace_(nullptr), fec_(nullptr), modeType_(mode), order_(order)
{
    if (!mesh_) { throw std::runtime_error("PhysicsWaveguideEigen: null Mesh"); }
    int dim = mesh_->Dimension();
    if (dim != 2) { throw std::runtime_error("PhysicsWaveguideEigen: only 2D meshes supported"); }

    fec_ = new mfem::H1_FECollection(order_, dim);
    fespace_ = new mfem::FiniteElementSpace(mesh_, fec_);
}
#endif

PhysicsWaveguideEigen::~PhysicsWaveguideEigen()
{
    delete stiffForm_;
    delete massForm_;
    delete fespace_;
    delete fec_;
}

void PhysicsWaveguideEigen::assembleStiffness()
{
#ifdef MFEM_USE_MPI
    stiffForm_ = new mfem::ParBilinearForm(fespace_);
#else
    stiffForm_ = new mfem::BilinearForm(fespace_);
#endif
    mfem::ConstantCoefficient one(1.0);
    stiffForm_->AddDomainIntegrator(new mfem::DiffusionIntegrator(one));

    // Assemble and, if TM with explicit boundaries, eliminate essential DOFs
    stiffForm_->Assemble();
    if (modeType_ == WaveguideModeType::TM)
    {
#ifdef MFEM_USE_MPI
        mfem::Array<int> essBdr;
        if (pmesh_->bdr_attributes.Size())
        {
            essBdr.SetSize(pmesh_->bdr_attributes.Max());
            essBdr = 0;
            pmesh_->MarkExternalBoundaries(essBdr);
            stiffForm_->EliminateEssentialBCDiag(essBdr, 1.0);
        }
#else
        mfem::Array<int> essBdr;
        if (mesh_->bdr_attributes.Size())
        {
            essBdr.SetSize(mesh_->bdr_attributes.Max());
            essBdr = 0;
            mesh_->MarkExternalBoundaries(essBdr);
            stiffForm_->EliminateEssentialBCDiag(essBdr, 1.0);
        }
#endif
    }
    stiffForm_->Finalize();
}

void PhysicsWaveguideEigen::assembleMass()
{
#ifdef MFEM_USE_MPI
    massForm_ = new mfem::ParBilinearForm(fespace_);
#else
    massForm_ = new mfem::BilinearForm(fespace_);
#endif
    mfem::ConstantCoefficient one(1.0);
    massForm_->AddDomainIntegrator(new mfem::MassIntegrator(one));
    massForm_->Assemble();
    if (modeType_ == WaveguideModeType::TM)
    {
#ifdef MFEM_USE_MPI
        mfem::Array<int> essBdr;
        if (pmesh_->bdr_attributes.Size())
        {
            essBdr.SetSize(pmesh_->bdr_attributes.Max());
            essBdr = 0;
            pmesh_->MarkExternalBoundaries(essBdr);
            massForm_->EliminateEssentialBCDiag(essBdr, std::numeric_limits<double>::min());
        }
#else
        mfem::Array<int> essBdr;
        if (mesh_->bdr_attributes.Size())
        {
            essBdr.SetSize(mesh_->bdr_attributes.Max());
            essBdr = 0;
            mesh_->MarkExternalBoundaries(essBdr);
            massForm_->EliminateEssentialBCDiag(essBdr, std::numeric_limits<double>::min());
        }
#endif
    }
    massForm_->Finalize();
}

std::vector<double> PhysicsWaveguideEigen::solveEigenvalues(int numModes)
{
    if (numModes <= 0) { throw std::runtime_error("PhysicsWaveguideEigen: numModes must be positive"); }

    eigenvalues_.clear();
    eigenvectors_.clear();

    assembleStiffness();
    assembleMass();

#ifdef MFEM_USE_MPI
    // Parallel path: form system matrices with essential BCs if TM
    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;

    if (modeType_ == WaveguideModeType::TM)
    {
        if (pmesh_->bdr_attributes.Size())
        {
            ess_bdr.SetSize(pmesh_->bdr_attributes.Max());
            ess_bdr = 0;
            pmesh_->MarkExternalBoundaries(ess_bdr);
            fespace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
        }
    }

    // Get parallel assembled operators (they already had essential DOFs handled above)
    mfem::HypreParMatrix *A = stiffForm_->ParallelAssemble();
    mfem::HypreParMatrix *M = massForm_->ParallelAssemble();

    // If TE (natural BC) the stiffness A is singular (constant nullspace).
    // Regularize by adding a rank-1 M-orthonormal projector: P = alpha * v v^T
    // where v is the constant true-dof vector normalized so that v^T M v = 1.
    if (modeType_ == WaveguideModeType::TE)
    {
        double alpha = 1e-8; // projector strength

        // Build a constant ParGridFunction and get its true-dofs vector
        mfem::ParGridFunction ones(fespace_);
        ones = 1.0;
        mfem::HypreParVector *tv = ones.GetTrueDofs();

        // Compute v^T M v and normalize v
        mfem::HypreParVector tmp = tv->CreateCompatibleVector();
        M->Mult(*tv, tmp);
        double mvv = mfem::InnerProduct(tv, &tmp);

        if (mvv <= 0.0 || !std::isfinite(mvv))
        {
            int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0)
            {
                std::cerr << "ERROR: Invalid M-norm for constant vector: mvv=" << mvv << std::endl;
            }
            delete tv;
            throw std::runtime_error("Cannot normalize constant vector for TE projector");
        }

        double scale = 1.0 / sqrt(mvv);
        *tv *= scale; // Now tv is M-orthonormal: tv^T M tv = 1

        int rank = 0, mpisize = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &mpisize);

        // Build the distributed rank-1 projector P = alpha * tv * tv^T
        // For HypreParMatrix, we need to construct diag and offd blocks.

        int local_size = tv->Size();
        HYPRE_BigInt *row_starts = fespace_->GetTrueDofOffsets();
        HYPRE_BigInt global_size = fespace_->GlobalTrueVSize();

        // Get local portion of the normalized vector
        mfem::Vector v_local(local_size);
        for (int i = 0; i < local_size; ++i)
        {
            v_local(i) = (*tv)(i);
        }

        // Gather all vector values across all ranks
        std::vector<double> v_global(global_size);
        std::vector<int> recvcounts(mpisize);
        std::vector<int> displs(mpisize);
        
        // Gather local sizes from all ranks
        MPI_Allgather(&local_size, 1, MPI_INT,
                      recvcounts.data(), 1, MPI_INT,
                      MPI_COMM_WORLD);

        // Compute displacements
        displs[0] = 0;
        for (int r = 1; r < mpisize; ++r)
        {
            displs[r] = displs[r-1] + recvcounts[r-1];
        }

        // Each rank sends its own local_size elements
        MPI_Allgatherv(v_local.GetData(), local_size, MPI_DOUBLE,
                       v_global.data(), recvcounts.data(), displs.data(),
                       MPI_DOUBLE, MPI_COMM_WORLD);

        // Determine the first column index for this rank
        HYPRE_BigInt first_col = displs[rank];

        // Construct diagonal block: P_diag(i,j) = alpha * v_local(i) * v_local(j)
        int *diag_i = new int[local_size + 1];
        int *diag_j = new int[local_size * local_size];
        double *diag_data = new double[local_size * local_size];
        
        diag_i[0] = 0;
        for (int i = 0; i < local_size; ++i)
        {
            diag_i[i+1] = diag_i[i] + local_size;
            for (int j = 0; j < local_size; ++j)
            {
                diag_j[diag_i[i] + j] = j;
                diag_data[diag_i[i] + j] = alpha * v_local(i) * v_local(j);
            }
        }

        // Determine off-diagonal columns (all columns NOT in local range)
        int offd_ncols = (int)(global_size - local_size);
        HYPRE_BigInt *col_map_offd = nullptr;
        int *offd_i = nullptr;
        int *offd_j = nullptr;
        double *offd_data = nullptr;

        if (offd_ncols > 0)
        {
            col_map_offd = new HYPRE_BigInt[offd_ncols];
            int col_idx = 0;
            for (HYPRE_BigInt j = 0; j < global_size; ++j)
            {
                if (j < first_col || j >= first_col + local_size)
                {
                    col_map_offd[col_idx++] = j;
                }
            }

            // Build off-diagonal CSR arrays
            offd_i = new int[local_size + 1];
            offd_j = new int[local_size * offd_ncols];
            offd_data = new double[local_size * offd_ncols];

            offd_i[0] = 0;
            for (int i = 0; i < local_size; ++i)
            {
                offd_i[i+1] = offd_i[i] + offd_ncols;
                int idx = 0;
                for (HYPRE_BigInt j = 0; j < global_size; ++j)
                {
                    if (j < first_col || j >= first_col + local_size)
                    {
                        offd_j[offd_i[i] + idx] = idx;
                        offd_data[offd_i[i] + idx] = alpha * v_local(i) * v_global[j];
                        idx++;
                    }
                }
            }
        }
        else
        {
            // Single rank case: no off-diagonal block
            offd_i = new int[local_size + 1];
            for (int i = 0; i <= local_size; ++i)
                offd_i[i] = 0;
        }

        // Create SparseMatrix wrappers (ownership transferred to SparseMatrix)
        mfem::SparseMatrix diag_mat(diag_i, diag_j, diag_data, local_size, local_size, 
                                     true, true, true);
        mfem::SparseMatrix offd_mat(offd_i, offd_j, offd_data, local_size, offd_ncols, 
                                     true, true, true);

        // Create HypreParMatrix from diag, offd, and col_map_offd
        // The HypreParMatrix constructor will copy data from SparseMatrix objects
        // and take ownership of col_map_offd
        mfem::HypreParMatrix *Ppar = new mfem::HypreParMatrix(MPI_COMM_WORLD,
                                                              global_size, global_size,
                                                              row_starts, row_starts,
                                                              &diag_mat, &offd_mat,
                                                              col_map_offd);

        // Add projector to A
        mfem::HypreParMatrix *Areg = mfem::Add(1.0, *A, 1.0, *Ppar);
        delete A;
        A = Areg;

        delete Ppar;
        delete tv;
    }

    // Preconditioner: HypreBoomerAMG for Laplacian
    mfem::HypreBoomerAMG *amg = new mfem::HypreBoomerAMG(*A);
    amg->SetPrintLevel(0);

    mfem::HypreLOBPCG *lobpcg = new mfem::HypreLOBPCG(MPI_COMM_WORLD);
    lobpcg->SetNumModes(numModes);
    lobpcg->SetRandomSeed(0);
    lobpcg->SetPreconditioner(*amg);
    lobpcg->SetMaxIter(200);
    lobpcg->SetTol(1e-8);
    lobpcg->SetPrecondUsageMode(1);
    lobpcg->SetPrintLevel(0);
    lobpcg->SetMassMatrix(*M);
    lobpcg->SetOperator(*A);

    lobpcg->Solve();

    mfem::Array<double> evals;
    lobpcg->GetEigenvalues(evals);

    for (int i = 0; i < evals.Size(); ++i)
    {
        eigenvalues_.push_back(evals[i]);

        // Extract eigenvector into a ParGridFunction
        mfem::ParGridFunction gf(fespace_);
        gf = lobpcg->GetEigenvector(i);

        // store local dofs for getEigenvectors()
        mfem::Vector local; local = gf;
        eigenvectors_.push_back(local);
    }

    delete lobpcg; delete amg; delete A; delete M;

#else
    // Serial path: prefer SLEPc (requires MFEM built with SLEPc). Use SLEPc wrapper.
#ifdef MFEM_USE_SLEPC
    mfem::Array<int> ess_bdr;
    mfem::Array<int> ess_tdof_list;
    if (modeType_ == WaveguideModeType::TM)
    {
        if (mesh_->bdr_attributes.Size())
        {
            ess_bdr.SetSize(mesh_->bdr_attributes.Max());
            ess_bdr = 0;
            mesh_->MarkExternalBoundaries(ess_bdr);
            fespace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
            std::cout << "TM (serial): essential true dofs = " << ess_tdof_list.Size() << std::endl;
        }
    }

    mfem::SparseMatrix *A = nullptr;
    mfem::SparseMatrix *M = nullptr;
    mfem::SparseMatrix *A_sys = nullptr;
    mfem::BilinearForm *sf = stiffForm_;
    mfem::BilinearForm *mf = massForm_;

    // FormSystemMatrix is available for serial BilinearForm too
    sf->FormSystemMatrix(ess_tdof_list, A);
    mf->FormSystemMatrix(ess_tdof_list, M);

    // If TE (natural BC) regularize serial A by adding rank-1 M-orthonormal projector
    if (modeType_ == WaveguideModeType::TE)
    {
        double alpha = 1e-8;
        int n = A->Height();
        mfem::Vector v(n);
        v = 1.0; // constant true-dof vector

        mfem::Vector tmp(n);
        M->Mult(v, tmp);
        double mvv = v * tmp; // v^T M v

        if (mvv > 0.0 && std::isfinite(mvv))
        {
            v *= 1.0 / sqrt(mvv); // now v^T M v == 1

            // Build rank-1 SparseMatrix P
            mfem::SparseMatrix *P = new mfem::SparseMatrix(n, n);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    double val = alpha * v(i) * v(j);
                    if (val != 0.0)
                    {
                        P->Add(i, j, val);
                    }
                }
            }
            P->Finalize();

            mfem::SparseMatrix *Areg = mfem::Add(1.0, *A, 1.0, *P);
            delete P;
            delete A;
            A = Areg;
        }
    }

    // Use MFEM SlepcEigenSolver wrapper
    mfem::SlepcEigenSolver slepc;
    slepc.SetNumModes(numModes);
    slepc.SetTolerances(1e-8, 200);
    // Request smallest magnitude eigenpairs (needed for Neumann/TE case)
    slepc.SetWhichEigenpairs(mfem::SlepcEigenSolver::SMALLEST_MAGNITUDE);
    slepc.SetOperator(*A, *M);
    slepc.Solve();

    mfem::Vector evals; slepc.GetEigenvalues(evals);
    for (int i = 0; i < evals.Size(); ++i)
    {
        eigenvalues_.push_back(evals(i));
        mfem::GridFunction gf(fespace_);
        gf = slepc.GetEigenvector(i);
        mfem::Vector local = gf;
        eigenvectors_.push_back(local);
    }

    delete A; delete M;
#else
    throw std::runtime_error("PhysicsWaveguideEigen: Serial eigen solver requires MFEM built with SLEPc (MFEM_USE_SLEPC). Please enable SLEPc or run in MPI mode with Hypre.");
#endif
#endif

    isSolved_ = true;
    return eigenvalues_;
}

std::vector<mfem::Vector> PhysicsWaveguideEigen::getEigenvectors() const
{
    if (!isSolved_) { throw std::runtime_error("PhysicsWaveguideEigen: call solveEigenvalues() first"); }
    return eigenvectors_;
}

#ifdef MFEM_USE_MPI
mfem::ParGridFunction* PhysicsWaveguideEigen::getModeShape(int modeIndex)
{
    if (!isSolved_) { throw std::runtime_error("PhysicsWaveguideEigen: call solveEigenvalues() first"); }
    if (modeIndex < 0 || modeIndex >= static_cast<int>(eigenvectors_.size()))
    { throw std::runtime_error("PhysicsWaveguideEigen: modeIndex out of range"); }

    mfem::ParGridFunction* gf = new mfem::ParGridFunction(fespace_);
    *gf = eigenvectors_[modeIndex];
    return gf;
}
#else
mfem::GridFunction* PhysicsWaveguideEigen::getModeShape(int modeIndex)
{
    if (!isSolved_) { throw std::runtime_error("PhysicsWaveguideEigen: call solveEigenvalues() first"); }
    if (modeIndex < 0 || modeIndex >= static_cast<int>(eigenvectors_.size()))
    { throw std::runtime_error("PhysicsWaveguideEigen: modeIndex out of range"); }

    mfem::GridFunction* gf = new mfem::GridFunction(fespace_);
    *gf = eigenvectors_[modeIndex];
    return gf;
}
#endif

} // namespace hpcfem
