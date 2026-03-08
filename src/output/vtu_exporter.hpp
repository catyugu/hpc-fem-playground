#ifndef MPFEM_VTU_EXPORTER_HPP
#define MPFEM_VTU_EXPORTER_HPP

#include "simulation_data.hpp"

#include <string>

namespace mpfem {

/**
 * @brief Exports sampled fields to VTU unstructured grid format.
 */
class VtuExporter {
public:
    /**
     * @brief Writes point-vertex VTU with scalar fields.
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
