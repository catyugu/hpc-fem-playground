#ifndef MPFEM_COMSOL_RESULT_READER_HPP
#define MPFEM_COMSOL_RESULT_READER_HPP

#include "simulation_data.hpp"

#include <vector>
#include <string>

namespace mpfem {

/**
 * @brief One row from COMSOL-style result text file.
 */
struct ComsolResultRow {
    Coordinate3D coordinate;
    double electricPotential = 0.0;
    double temperature = 0.0;
    double displacement = 0.0;
};

/**
 * @brief Reads COMSOL text result format used by busbar reference.
 */
class ComsolResultReader {
public:
    /**
     * @brief Parses rows from result text file.
     * @param filePath Input result path.
     * @param rows Parsed rows.
     * @param errorMessage Error details on failure.
     * @return True when parsing succeeds.
     */
    static bool readFromFile(const std::string &filePath,
                             std::vector<ComsolResultRow> &rows,
                             std::string &errorMessage);
};

} // namespace mpfem

#endif
