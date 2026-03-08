#ifndef MPFEM_COMSOL_RESULT_READER_HPP
#define MPFEM_COMSOL_RESULT_READER_HPP

#include "simulation_data.hpp"

#include <string>
#include <vector>

namespace mpfem {

/**
 * @brief Reads COMSOL reference results from text file.
 */
class ComsolResultReader {
public:
    /**
     * @brief Parses COMSOL result file.
     */
    static void readFromFile(const std::string &filePath, std::vector<ComsolResultRow> &rows);
};

} // namespace mpfem

#endif