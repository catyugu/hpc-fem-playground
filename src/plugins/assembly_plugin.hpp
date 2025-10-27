#ifndef MFEM_ASSEMBLY_PLUGIN_HPP
#define MFEM_ASSEMBLY_PLUGIN_HPP

#include "plugin_base.hpp"

namespace mfem_plugins
{

    // Assembly strategy enumeration
    enum class AssemblyStrategy
    {
        FULL,              // Standard full assembly
        PARTIAL,           // Partial assembly (element matrices stored)
        MATRIX_FREE,       // Matrix-free (on-the-fly computation)
        ELEMENT_BY_ELEMENT // Element-by-element with custom storage
    };

    // Base class for assembly plugins
    class AssemblyPlugin : public PluginBase
    {
    protected:
        AssemblyStrategy strategy_;
        mfem::FiniteElementSpace *fes_;

    public:
        AssemblyPlugin(const std::string &name, AssemblyStrategy strategy)
            : PluginBase(name), strategy_(strategy), fes_(nullptr) {}

        virtual ~AssemblyPlugin() {}

        AssemblyStrategy strategy() const { return strategy_; }

        virtual void set_fespace(mfem::FiniteElementSpace *fes) { fes_ = fes; }

        // Assemble bilinear form (stiffness/mass matrix)
        virtual void assemble(mfem::BilinearForm &form) = 0;

        // Apply operator: y = A*x
        virtual void mult(const mfem::Vector &x, mfem::Vector &y) = 0;

        // Get sparse matrix (if applicable)
        virtual mfem::SparseMatrix *get_matrix() { return nullptr; }

        // Memory usage estimation
        virtual size_t memory_usage() const = 0;
    };

    // Standard full assembly plugin
    class FullAssemblyPlugin : public AssemblyPlugin
    {
    private:
        mfem::SparseMatrix *matrix_;

    public:
        FullAssemblyPlugin()
            : AssemblyPlugin("Full Assembly", AssemblyStrategy::FULL),
              matrix_(nullptr) {}

        virtual ~FullAssemblyPlugin() {}

        virtual void assemble(mfem::BilinearForm &form) override
        {
            Timer timer;
            timer.start();

            form.Assemble();
            form.Finalize();
            matrix_ = &form.SpMat();

            metrics_.setup_time_ms = timer.stop();
            metrics_.memory_bytes = memory_usage();
        }

        virtual void mult(const mfem::Vector &x, mfem::Vector &y) override
        {
            if (matrix_)
            {
                matrix_->Mult(x, y);
            }
        }

        virtual mfem::SparseMatrix *get_matrix() override
        {
            return matrix_;
        }

        virtual size_t memory_usage() const override
        {
            if (!matrix_)
                return 0;
            // Approximate: nnz * (sizeof(double) + sizeof(int)) + row_ptr
            return matrix_->NumNonZeroElems() * (sizeof(double) + sizeof(int)) +
                   (matrix_->Height() + 1) * sizeof(int);
        }
    };

    // Partial assembly plugin
    class PartialAssemblyPlugin : public AssemblyPlugin
    {
    private:
        mfem::BilinearForm *form_;

    public:
        PartialAssemblyPlugin()
            : AssemblyPlugin("Partial Assembly", AssemblyStrategy::PARTIAL),
              form_(nullptr) {}

        virtual ~PartialAssemblyPlugin() {}

        virtual void assemble(mfem::BilinearForm &form) override
        {
            Timer timer;
            timer.start();

            form.SetAssemblyLevel(mfem::AssemblyLevel::PARTIAL);
            form.Assemble();
            form_ = &form;

            metrics_.setup_time_ms = timer.stop();
            metrics_.memory_bytes = memory_usage();
        }

        virtual void mult(const mfem::Vector &x, mfem::Vector &y) override
        {
            if (form_)
            {
                form_->Mult(x, y);
            }
        }

        virtual size_t memory_usage() const override
        {
            if (!form_ || !fes_)
                return 0;
            // Estimate: number of elements * element matrix size * sizeof(double)
            int ne = fes_->GetNE();
            int dofs_per_elem = fes_->GetFE(0)->GetDof();
            return ne * dofs_per_elem * dofs_per_elem * sizeof(double);
        }
    };

    // Matrix-free assembly plugin
    class MatrixFreePlugin : public AssemblyPlugin
    {
    private:
        mfem::BilinearForm *form_;

    public:
        MatrixFreePlugin()
            : AssemblyPlugin("Matrix-Free", AssemblyStrategy::MATRIX_FREE),
              form_(nullptr) {}

        virtual ~MatrixFreePlugin() {}

        virtual void assemble(mfem::BilinearForm &form) override
        {
            Timer timer;
            timer.start();

            // Matrix-free: no assembly, just store form reference
            form.SetAssemblyLevel(mfem::AssemblyLevel::NONE);
            form.Assemble();
            form_ = &form;

            metrics_.setup_time_ms = timer.stop();
            metrics_.memory_bytes = memory_usage();
        }

        virtual void mult(const mfem::Vector &x, mfem::Vector &y) override
        {
            if (form_)
            {
                form_->Mult(x, y);
            }
        }

        virtual size_t memory_usage() const override
        {
            // Minimal memory: just coefficient storage
            return 1024; // Placeholder estimate
        }
    };

    // Element-by-element assembly with custom storage
    class ElementByElementPlugin : public AssemblyPlugin
    {
    private:
        std::vector<mfem::DenseMatrix> elem_matrices_;
        mfem::BilinearForm *form_;

    public:
        ElementByElementPlugin()
            : AssemblyPlugin("Element-by-Element", AssemblyStrategy::ELEMENT_BY_ELEMENT),
              form_(nullptr) {}

        virtual ~ElementByElementPlugin() {}

        virtual void assemble(mfem::BilinearForm &form) override
        {
            Timer timer;
            timer.start();

            form_ = &form;
            if (!fes_)
                return;

            int ne = fes_->GetNE();
            elem_matrices_.reserve(ne);

            // Assemble element matrices
            mfem::Array<int> vdofs;
            mfem::ElementTransformation *eltrans;

            for (int i = 0; i < ne; i++)
            {
                const mfem::FiniteElement &fe = *fes_->GetFE(i);
                fes_->GetElementVDofs(i, vdofs);

                mfem::DenseMatrix elmat;
                elmat.SetSize(vdofs.Size());
                elmat = 0.0;

                eltrans = fes_->GetElementTransformation(i);

                // Get integrators from form and assemble element matrix
                // Note: This is simplified - real implementation would iterate integrators
                elem_matrices_.push_back(elmat);
            }

            metrics_.setup_time_ms = timer.stop();
            metrics_.memory_bytes = memory_usage();
        }

        virtual void mult(const mfem::Vector &x, mfem::Vector &y) override
        {
            if (!fes_ || elem_matrices_.empty())
                return;

            y = 0.0;
            mfem::Array<int> vdofs;
            mfem::Vector x_elem, y_elem;

            for (size_t i = 0; i < elem_matrices_.size(); i++)
            {
                fes_->GetElementVDofs(i, vdofs);
                x.GetSubVector(vdofs, x_elem);

                y_elem.SetSize(vdofs.Size());
                elem_matrices_[i].Mult(x_elem, y_elem);

                y.AddElementVector(vdofs, y_elem);
            }
        }

        virtual size_t memory_usage() const override
        {
            size_t total = 0;
            for (const auto &mat : elem_matrices_)
            {
                total += mat.Height() * mat.Width() * sizeof(double);
            }
            return total;
        }
    };

    // Factory function
    inline std::unique_ptr<AssemblyPlugin> create_assembly_plugin(AssemblyStrategy strategy)
    {
        switch (strategy)
        {
        case AssemblyStrategy::FULL:
            return std::unique_ptr<AssemblyPlugin>(new FullAssemblyPlugin());
        case AssemblyStrategy::PARTIAL:
            return std::unique_ptr<AssemblyPlugin>(new PartialAssemblyPlugin());
        case AssemblyStrategy::MATRIX_FREE:
            return std::unique_ptr<AssemblyPlugin>(new MatrixFreePlugin());
        case AssemblyStrategy::ELEMENT_BY_ELEMENT:
            return std::unique_ptr<AssemblyPlugin>(new ElementByElementPlugin());
        default:
            return std::unique_ptr<AssemblyPlugin>(new FullAssemblyPlugin());
        }
    }

} // namespace mfem_plugins

#endif // MFEM_ASSEMBLY_PLUGIN_HPP
