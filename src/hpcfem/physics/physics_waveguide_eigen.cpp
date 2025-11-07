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
    // print mesh extents for diagnostics
    {
        mfem::Vector pmin, pmax;
        pmesh_->GetBoundingBox(pmin, pmax);
        int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0)
        {
            std::cout << "Mesh bbox: [" << pmin(0) << "," << pmin(1) << "] - ["
                      << pmax(0) << "," << pmax(1) << "]" << std::endl;
        }
    }

    if (modeType_ == WaveguideModeType::TM)
    {
        if (pmesh_->bdr_attributes.Size())
        {
            ess_bdr.SetSize(pmesh_->bdr_attributes.Max());
            ess_bdr = 0;
            pmesh_->MarkExternalBoundaries(ess_bdr);
            fespace_->GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
            int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0)
            {
                std::cout << "TM: essential true dofs = " << ess_tdof_list.Size() << std::endl;
            }
        }
    }

    // Get parallel assembled operators (they already had essential DOFs handled above)
    mfem::HypreParMatrix *A = stiffForm_->ParallelAssemble();
    mfem::HypreParMatrix *M = massForm_->ParallelAssemble();

    // Diagnostic: print diagonal stats for A and M to check scaling
    {
        mfem::Vector A_diag, M_diag;
        A->GetDiag(A_diag);
        M->GetDiag(M_diag);
        double A_min=1e300, A_max=-1e300, A_sum=0.0;
        double M_min=1e300, M_max=-1e300, M_sum=0.0;
        for (int i=0;i<A_diag.Size();++i){ double v=A_diag(i); A_min=std::min(A_min,v); A_max=std::max(A_max,v); A_sum+=v; }
        for (int i=0;i<M_diag.Size();++i){ double v=M_diag(i); M_min=std::min(M_min,v); M_max=std::max(M_max,v); M_sum+=v; }
        double A_mean = A_sum / A_diag.Size();
        double M_mean = M_sum / M_diag.Size();
        int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rank == 0)
        {
            std::cout << "A diag: min="<<A_min<<" max="<<A_max<<" mean="<<A_mean<<"\n";
            std::cout << "M diag: min="<<M_min<<" max="<<M_max<<" mean="<<M_mean<<"\n";
        }
    }

    // If TE (natural BC) the stiffness A is singular (constant nullspace).
    // Prefer a rank-1 M-orthonormal projector regularization (option 1).
    if (modeType_ == WaveguideModeType::TE)
    {
        double alpha = 1e-8; // base strength for projector

        // Build a constant ParGridFunction and get its true-dofs vector
        mfem::ParGridFunction ones(fespace_);
        ones = 1.0;
        mfem::HypreParVector *tv = ones.GetTrueDofs();

        // tmp = M * tv
        mfem::HypreParVector tmp = tv->CreateCompatibleVector();
        M->Mult(*tv, tmp);
        double mvv = mfem::InnerProduct(tv, &tmp);

        int mpisize = 1; MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
        if (mvv <= 0.0 || !std::isfinite(mvv) || mpisize != 1)
        {
            // Fallback: simple Tikhonov shift with the mass matrix.
            mfem::HypreParMatrix *Areg = mfem::Add(1.0, *A, alpha, *M);
            delete A;
            A = Areg;
        }
        else
        {
            // Single-rank MPI: construct a full rank-1 CSR on rank 0 and
            // convert it to a HypreParMatrix via the rank-0 CSR constructor.
            std::cout << "(parallel single-rank) TE: building rank-1 M-orthonormal projector, mvv=" << mvv << "\n";
            double scale = 1.0 / sqrt(mvv);

            // Gather global vector on rank 0
            mfem::Vector *gvp = tv->GlobalVector(); // full vector on rank 0
            mfem::Vector gv = *gvp;
            delete gvp;
            int n = gv.Size();

            // Build CSR arrays (dense outer-product) on rank 0
            int *I = new int[n+1];
            HYPRE_BigInt *Jbig = new HYPRE_BigInt[n*n];
            double *data = new double[n*n];
            I[0] = 0;
            for (int i = 0; i < n; ++i)
            {
                I[i+1] = I[i] + n;
                for (int j = 0; j < n; ++j)
                {
                    Jbig[I[i] + j] = (HYPRE_BigInt) j;
                    data[I[i] + j] = alpha * (gv(i)*scale) * (gv(j)*scale);
                }
            }

            // Create SparseMatrix from CSR arrays (ownership transferred)
            // Create HypreParMatrix from the CSR arrays (works with assumed
            // partition mode). Ownership of I/Jbig/data is transferred to
            // the HypreParMatrix constructor, so we don't free them here.
            HYPRE_BigInt *row_starts = fespace_->GetTrueDofOffsets();
            HYPRE_BigInt *col_starts = row_starts;
            mfem::HypreParMatrix *Ppar = new mfem::HypreParMatrix(MPI_COMM_WORLD,
                                                                 n, /* local rows */
                                                                 (HYPRE_BigInt) n, (HYPRE_BigInt) n,
                                                                 I, Jbig, data,
                                                                 row_starts, col_starts);

            // Add projector to A
            mfem::HypreParMatrix *Areg = mfem::Add(1.0, *A, 1.0, *Ppar);
            delete A;
            A = Areg;

            delete Ppar;
        }
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

    int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        std::cout << "Computed " << evals.Size() << " eigenvalues (k_c^2)" << std::endl;
        for (int ii = 0; ii < evals.Size(); ++ii)
        {
            std::cout << "  Mode " << ii << ": k_c^2 = " << evals[ii] << std::endl;
        }
    }

    for (int i = 0; i < evals.Size(); ++i)
    {
        eigenvalues_.push_back(evals[i]);

        // Extract eigenvector into a ParGridFunction
        mfem::ParGridFunction gf(fespace_);
        gf = lobpcg->GetEigenvector(i);

        // Diagnostic: compute Rayleigh quotient r = (v^T A v) / (v^T M v)
        {
            mfem::HypreParVector *tv = gf.GetTrueDofs();
            mfem::HypreParVector tmpA(*A), tmpM(*M);
            A->Mult(*tv, tmpA);
            M->Mult(*tv, tmpM);
            double num = mfem::InnerProduct(tv, &tmpA);
            double den = mfem::InnerProduct(tv, &tmpM);
            int rank = 0; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            if (rank == 0)
            {
                std::cout << "Mode " << i << ": eigenvalue=" << evals[i]
                          << " Rayleigh=" << (num/den) << std::endl;
            }
            delete tv;
        }

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

    // Diagnostic: print diagonal stats for A and M (serial)
    {
        mfem::Vector A_diag, M_diag;
        A->GetDiag(A_diag);
        M->GetDiag(M_diag);
        double A_min=1e300, A_max=-1e300, A_sum=0.0;
        double M_min=1e300, M_max=-1e300, M_sum=0.0;
        for (int i=0;i<A_diag.Size();++i){ double v=A_diag(i); A_min=std::min(A_min,v); A_max=std::max(A_max,v); A_sum+=v; }
        for (int i=0;i<M_diag.Size();++i){ double v=M_diag(i); M_min=std::min(M_min,v); M_max=std::max(M_max,v); M_sum+=v; }
        double A_mean = A_sum / A_diag.Size();
        double M_mean = M_sum / M_diag.Size();
        std::cout << "(serial) A diag: min="<<A_min<<" max="<<A_max<<" mean="<<A_mean<<"\n";
        std::cout << "(serial) M diag: min="<<M_min<<" max="<<M_max<<" mean="<<M_mean<<"\n";
    }

    // If TE (natural BC) regularize serial A by adding a small multiple of M
    if (modeType_ == WaveguideModeType::TE)
    {
        // Implement option (1): add a rank-1 M-orthonormal projector
        // P = alpha * v * v^T where v is the constant vector normalized
        // in the M-inner-product: v^T M v = 1.
        double alpha = 1e-8;
        int n = A->Height();
        mfem::Vector v(n);
        v = 1.0; // constant true-dof vector

        // tmp = M * v
        mfem::Vector tmp(n);
        M->Mult(v, tmp);
        double mvv = v * tmp; // v^T M v

        if (mvv <= 0.0 || !std::isfinite(mvv))
        {
            // fallback to diagonal regularization if something unexpected
            mfem::SparseMatrix *Areg = mfem::Add(1.0, *A, alpha, *M);
            delete A;
            A = Areg;
        }
        else
        {
                std::cout << "(serial) TE: building rank-1 M-orthonormal projector, mvv=" << mvv << "\n";
            v *= 1.0 / sqrt(mvv); // now v^T M v == 1

            // Build a (dense-in-sparsity) rank-1 SparseMatrix P locally
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
    std::cout << "Computed " << evals.Size() << " eigenvalues (k_c^2)" << std::endl;
    for (int i = 0; i < evals.Size(); ++i)
    {
        eigenvalues_.push_back(evals(i));
        mfem::GridFunction gf(fespace_);
        gf = slepc.GetEigenvector(i);
        mfem::Vector local = gf;
        // Diagnostic: compute Rayleigh quotient r = (v^T A v) / (v^T M v)
        {
            mfem::Vector tv; gf.GetTrueDofs(tv);
            mfem::Vector tmpA(tv.Size()), tmpM(tv.Size());
            A->Mult(tv, tmpA);
            M->Mult(tv, tmpM);
            double num = tv * tmpA;
            double den = tv * tmpM;
            std::cout << "(serial) Mode " << i << ": eigenvalue=" << evals(i)
                      << " Rayleigh=" << (num/den) << std::endl;
        }
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
