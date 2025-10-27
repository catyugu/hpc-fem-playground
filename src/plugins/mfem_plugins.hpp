#ifndef MFEM_PLUGINS_HPP
#define MFEM_PLUGINS_HPP

// Main header file for MFEM optimization plugins
// Include this file to access all plugin functionality

#include "plugin_base.hpp"
#include "assembly_plugin.hpp"
#include "solver_plugin.hpp"

namespace mfem_plugins
{

    // Version information
    constexpr const char *PLUGIN_VERSION = "1.0.0";
    constexpr const char *PLUGIN_NAME = "MFEM Performance Plugins";

    // Print plugin system information
    inline void print_plugin_info()
    {
        std::cout << "========================================" << std::endl;
        std::cout << PLUGIN_NAME << " v" << PLUGIN_VERSION << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Available Assembly Strategies:" << std::endl;
        std::cout << "  - FULL: Standard full assembly" << std::endl;
        std::cout << "  - PARTIAL: Partial assembly (element matrices)" << std::endl;
        std::cout << "  - MATRIX_FREE: Matrix-free operators" << std::endl;
        std::cout << "  - ELEMENT_BY_ELEMENT: Custom element storage" << std::endl;
        std::cout << "\nAvailable Solver Types:" << std::endl;
        std::cout << "  - CG: Conjugate Gradient" << std::endl;
        std::cout << "  - GMRES: Generalized Minimal Residual" << std::endl;
        std::cout << "  - DIRECT: Direct solver (UMFPACK)" << std::endl;
        std::cout << "\nAvailable Preconditioners:" << std::endl;
        std::cout << "  - NONE: No preconditioning" << std::endl;
        std::cout << "  - JACOBI: Jacobi (diagonal)" << std::endl;
        std::cout << "  - GAUSS_SEIDEL: Gauss-Seidel" << std::endl;
        std::cout << "  - AMG: Algebraic Multigrid (HYPRE)" << std::endl;
        std::cout << "========================================" << std::endl;
    }

} // namespace mfem_plugins

#endif // MFEM_PLUGINS_HPP
