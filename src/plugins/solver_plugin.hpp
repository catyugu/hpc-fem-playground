#ifndef MFEM_SOLVER_PLUGIN_HPP
#define MFEM_SOLVER_PLUGIN_HPP

#include "plugin_base.hpp"

namespace mfem_plugins {

// Solver type enumeration
enum class SolverType {
    CG,             // Conjugate Gradient
    GMRES,          // GMRES
    MINRES,         // MINRES
    BICGSTAB,       // BiCGSTAB
    DIRECT          // Direct solver (UMFPACK, etc.)
};

// Preconditioner type enumeration
enum class PreconditionerType {
    NONE,           // No preconditioning
    JACOBI,         // Jacobi (diagonal)
    GAUSS_SEIDEL,   // Gauss-Seidel
    AMG,            // Algebraic Multigrid (HYPRE BoomerAMG)
    ILU,            // Incomplete LU
    CUSTOM          // Custom preconditioner
};

// Base class for solver plugins
class SolverPlugin : public PluginBase {
protected:
    SolverType solver_type_;
    PreconditionerType precond_type_;
    double rel_tol_;
    double abs_tol_;
    int max_iter_;
    int print_level_;
    
    mfem::Solver* solver_;
    mfem::Solver* precond_;
    
public:
    SolverPlugin(const std::string& name, SolverType stype, PreconditionerType ptype)
        : PluginBase(name), solver_type_(stype), precond_type_(ptype),
          rel_tol_(1e-9), abs_tol_(1e-12), max_iter_(1000), print_level_(0),
          solver_(nullptr), precond_(nullptr) {}
    
    virtual ~SolverPlugin() {
        delete solver_;
        delete precond_;
    }
    
    void set_tolerance(double rel_tol, double abs_tol = 0.0) {
        rel_tol_ = rel_tol;
        abs_tol_ = abs_tol;
    }
    
    void set_max_iterations(int max_iter) { max_iter_ = max_iter; }
    void set_print_level(int level) { print_level_ = level; }
    
    // Setup solver with operator
    virtual void setup(mfem::Operator& op) = 0;
    
    // Solve: A*x = b
    virtual void solve(const mfem::Vector& b, mfem::Vector& x) = 0;
    
    SolverType solver_type() const { return solver_type_; }
    PreconditionerType precond_type() const { return precond_type_; }
};

// CG solver plugin
class CGSolverPlugin : public SolverPlugin {
private:
    mfem::CGSolver* cg_solver_;
    
public:
    CGSolverPlugin(PreconditionerType ptype = PreconditionerType::NONE)
        : SolverPlugin("CG Solver", SolverType::CG, ptype),
          cg_solver_(nullptr) {}
    
    virtual ~CGSolverPlugin() {}
    
    virtual void setup(mfem::Operator& op) override {
        Timer timer;
        timer.start();
        
        cg_solver_ = new mfem::CGSolver();
        cg_solver_->SetRelTol(rel_tol_);
        cg_solver_->SetAbsTol(abs_tol_);
        cg_solver_->SetMaxIter(max_iter_);
        cg_solver_->SetPrintLevel(print_level_);
        
        // Setup preconditioner
        if (precond_type_ == PreconditionerType::JACOBI) {
            precond_ = new mfem::DSmoother();
        } else if (precond_type_ == PreconditionerType::GAUSS_SEIDEL) {
            precond_ = new mfem::GSSmoother();
        }
        
        if (precond_) {
            precond_->SetOperator(op);
            cg_solver_->SetPreconditioner(*precond_);
        }
        
        cg_solver_->SetOperator(op);
        solver_ = cg_solver_;
        
        metrics_.setup_time_ms = timer.stop();
    }
    
    virtual void solve(const mfem::Vector& b, mfem::Vector& x) override {
        if (!cg_solver_) return;
        
        Timer timer;
        timer.start();
        
        cg_solver_->Mult(b, x);
        
        metrics_.solve_time_ms = timer.stop();
        metrics_.total_time_ms = metrics_.setup_time_ms + metrics_.solve_time_ms;
        metrics_.iterations = cg_solver_->GetNumIterations();
        metrics_.residual = cg_solver_->GetFinalNorm();
    }
};

// GMRES solver plugin
class GMRESSolverPlugin : public SolverPlugin {
private:
    mfem::GMRESSolver* gmres_solver_;
    int krylov_dim_;
    
public:
    GMRESSolverPlugin(PreconditionerType ptype = PreconditionerType::NONE, int kdim = 50)
        : SolverPlugin("GMRES Solver", SolverType::GMRES, ptype),
          gmres_solver_(nullptr), krylov_dim_(kdim) {}
    
    virtual ~GMRESSolverPlugin() {}
    
    void set_krylov_dim(int kdim) { krylov_dim_ = kdim; }
    
    virtual void setup(mfem::Operator& op) override {
        Timer timer;
        timer.start();
        
        gmres_solver_ = new mfem::GMRESSolver();
        gmres_solver_->SetRelTol(rel_tol_);
        gmres_solver_->SetAbsTol(abs_tol_);
        gmres_solver_->SetMaxIter(max_iter_);
        gmres_solver_->SetPrintLevel(print_level_);
        gmres_solver_->SetKDim(krylov_dim_);
        
        // Setup preconditioner
        if (precond_type_ == PreconditionerType::JACOBI) {
            precond_ = new mfem::DSmoother();
        } else if (precond_type_ == PreconditionerType::GAUSS_SEIDEL) {
            precond_ = new mfem::GSSmoother();
        }
        
        if (precond_) {
            precond_->SetOperator(op);
            gmres_solver_->SetPreconditioner(*precond_);
        }
        
        gmres_solver_->SetOperator(op);
        solver_ = gmres_solver_;
        
        metrics_.setup_time_ms = timer.stop();
    }
    
    virtual void solve(const mfem::Vector& b, mfem::Vector& x) override {
        if (!gmres_solver_) return;
        
        Timer timer;
        timer.start();
        
        gmres_solver_->Mult(b, x);
        
        metrics_.solve_time_ms = timer.stop();
        metrics_.total_time_ms = metrics_.setup_time_ms + metrics_.solve_time_ms;
        metrics_.iterations = gmres_solver_->GetNumIterations();
        metrics_.residual = gmres_solver_->GetFinalNorm();
    }
};

// AMG preconditioned CG solver plugin
class AMGCGSolverPlugin : public SolverPlugin {
private:
    mfem::CGSolver* cg_solver_;
    mfem::HypreBoomerAMG* amg_precond_;
    
public:
    AMGCGSolverPlugin()
        : SolverPlugin("AMG-CG Solver", SolverType::CG, PreconditionerType::AMG),
          cg_solver_(nullptr), amg_precond_(nullptr) {}
    
    virtual ~AMGCGSolverPlugin() {
        delete amg_precond_;
    }
    
    virtual void setup(mfem::Operator& op) override {
        Timer timer;
        timer.start();
        
        // Setup AMG preconditioner
        mfem::HypreParMatrix* hypre_op = dynamic_cast<mfem::HypreParMatrix*>(&op);
        if (hypre_op) {
            amg_precond_ = new mfem::HypreBoomerAMG(*hypre_op);
            amg_precond_->SetPrintLevel(0);
        } else {
            // Fallback for serial matrices
            mfem::SparseMatrix* sparse_op = dynamic_cast<mfem::SparseMatrix*>(&op);
            if (sparse_op) {
                // Use simpler preconditioner for serial case
                precond_ = new mfem::GSSmoother(*sparse_op);
            }
        }
        
        // Setup CG solver
        cg_solver_ = new mfem::CGSolver();
        cg_solver_->SetRelTol(rel_tol_);
        cg_solver_->SetAbsTol(abs_tol_);
        cg_solver_->SetMaxIter(max_iter_);
        cg_solver_->SetPrintLevel(print_level_);
        
        if (amg_precond_) {
            cg_solver_->SetPreconditioner(*amg_precond_);
        } else if (precond_) {
            cg_solver_->SetPreconditioner(*precond_);
        }
        
        cg_solver_->SetOperator(op);
        solver_ = cg_solver_;
        
        metrics_.setup_time_ms = timer.stop();
    }
    
    virtual void solve(const mfem::Vector& b, mfem::Vector& x) override {
        if (!cg_solver_) return;
        
        Timer timer;
        timer.start();
        
        cg_solver_->Mult(b, x);
        
        metrics_.solve_time_ms = timer.stop();
        metrics_.total_time_ms = metrics_.setup_time_ms + metrics_.solve_time_ms;
        metrics_.iterations = cg_solver_->GetNumIterations();
        metrics_.residual = cg_solver_->GetFinalNorm();
    }
};

// Direct solver plugin (UMFPACK)
class DirectSolverPlugin : public SolverPlugin {
private:
#ifdef MFEM_USE_SUITESPARSE
    mfem::UMFPackSolver* direct_solver_;
#endif
    
public:
    DirectSolverPlugin()
        : SolverPlugin("Direct Solver (UMFPACK)", SolverType::DIRECT, PreconditionerType::NONE)
#ifdef MFEM_USE_SUITESPARSE
          , direct_solver_(nullptr)
#endif
    {}
    
    virtual ~DirectSolverPlugin() {
#ifdef MFEM_USE_SUITESPARSE
        delete direct_solver_;
#endif
    }
    
    virtual void setup(mfem::Operator& op) override {
        Timer timer;
        timer.start();
        
#ifdef MFEM_USE_SUITESPARSE
        direct_solver_ = new mfem::UMFPackSolver();
        direct_solver_->SetOperator(op);
        solver_ = direct_solver_;
#else
        std::cerr << "UMFPACK not available. Compile MFEM with SuiteSparse." << std::endl;
#endif
        
        metrics_.setup_time_ms = timer.stop();
    }
    
    virtual void solve(const mfem::Vector& b, mfem::Vector& x) override {
#ifdef MFEM_USE_SUITESPARSE
        if (!direct_solver_) return;
        
        Timer timer;
        timer.start();
        
        direct_solver_->Mult(b, x);
        
        metrics_.solve_time_ms = timer.stop();
        metrics_.total_time_ms = metrics_.setup_time_ms + metrics_.solve_time_ms;
        metrics_.iterations = 1;  // Direct solve
        metrics_.residual = 0.0;  // Exact (within machine precision)
#endif
    }
};

// Factory function
inline std::unique_ptr<SolverPlugin> create_solver_plugin(
    SolverType stype, 
    PreconditionerType ptype = PreconditionerType::NONE) {
    
    switch (stype) {
        case SolverType::CG:
            if (ptype == PreconditionerType::AMG) {
                return std::unique_ptr<SolverPlugin>(new AMGCGSolverPlugin());
            }
            return std::unique_ptr<SolverPlugin>(new CGSolverPlugin(ptype));
            
        case SolverType::GMRES:
            return std::unique_ptr<SolverPlugin>(new GMRESSolverPlugin(ptype));
            
        case SolverType::DIRECT:
            return std::unique_ptr<SolverPlugin>(new DirectSolverPlugin());
            
        default:
            return std::unique_ptr<SolverPlugin>(new CGSolverPlugin(ptype));
    }
}

} // namespace mfem_plugins

#endif // MFEM_SOLVER_PLUGIN_HPP
