#include "comsol_result_reader.hpp"

#include <fstream>
#include <sstream>

namespace mpfem {

bool ComsolResultReader::readFromFile(const std::string &filePath,
                                      std::vector<ComsolResultRow> &rows,
                                      std::string &errorMessage)
{
    rows.clear();
    errorMessage.clear();

    std::ifstream input(filePath.c_str());
    if (!input.is_open()) {
        errorMessage = "Cannot open COMSOL result file: " + filePath;
        return false;
    }

    std::string line;
    while (std::getline(input, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '%') {
            continue;
        }
        if (line[0] == 'x' || line[0] == 'X') {
            continue;
        }

        std::stringstream stream(line);
        ComsolResultRow row;
        if (!(stream >> row.coordinate.x
                     >> row.coordinate.y
                     >> row.coordinate.z
                     >> row.electricPotential
                     >> row.temperature
                     >> row.displacement)) {
            continue;
        }

        rows.push_back(row);
    }

    if (rows.empty()) {
        errorMessage = "No numeric rows parsed from COMSOL result file: " + filePath;
        return false;
    }

    return true;
}

} // namespace mpfem
