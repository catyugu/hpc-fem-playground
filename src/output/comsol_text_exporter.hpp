#ifndef MPFEM_COMSOL_TEXT_EXPORTER_HPP
#define MPFEM_COMSOL_TEXT_EXPORTER_HPP

#include "simulation_data.hpp"

#include <string>

namespace mpfem {

/**
 * @brief Exports results in COMSOL text format.
 */
class ComsolTextExporter {
public:
    /**
     * @brief Writes results to file.
     */
    static void write(const std::string &filePath, const CoupledFieldResult &result);
};

} // namespace mpfem

#endif