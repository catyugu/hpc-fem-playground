#ifndef MPFEM_VTU_EXPORTER_HPP
#define MPFEM_VTU_EXPORTER_HPP

#include "simulation_data.hpp"

#include <string>

namespace mpfem {

/**
 * @brief Exports results in VTU format.
 */
class VtuExporter {
public:
    /**
     * @brief Writes results to file.
     */
    static void write(const std::string &filePath, const CoupledFieldResult &result);
};

} // namespace mpfem

#endif