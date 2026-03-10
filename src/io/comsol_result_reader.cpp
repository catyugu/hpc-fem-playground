#include "comsol_result_reader.hpp"
#include "logger.hpp"

#include <fstream>
#include <sstream>

namespace mpfem {

void ComsolResultReader::readFromFile(const std::string &filePath, std::vector<ComsolResultRow> &rows)
{
    rows.clear();

    std::ifstream file(filePath);
    Check(file.is_open(), "Cannot open COMSOL result file: " + filePath);

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '%') {
            continue;
        }

        std::stringstream ss(line);
        ComsolResultRow row;
        ss >> row.x >> row.y >> row.z >> row.potential >> row.temperature >> row.displacement;
        rows.push_back(row);
    }

    Check(!rows.empty(), "No numeric rows parsed from COMSOL result file: " + filePath);
}

} // namespace mpfem