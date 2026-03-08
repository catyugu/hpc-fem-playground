#ifndef MPFEM_COMSOL_TEXT_EXPORTER_HPP
#define MPFEM_COMSOL_TEXT_EXPORTER_HPP

#include "simulation_data.hpp"

#include <string>

namespace mpfem {

/**
 * @brief Exports simulation data in COMSOL-like text format.
 */
class ComsolTextExporter {
public:
    /**
     * @brief Writes `x y z V T disp` table with metadata header.
     * @param outputPath Target file path.
     * @param result Coupled field result.
     * @param errorMessage Error details on failure.
     * @return True when write succeeds.
     */
    static bool write(const std::string &outputPath,
                      const CoupledFieldResult &result,
                      std::string &errorMessage);
};

} // namespace mpfem

#endif
