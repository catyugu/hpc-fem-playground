#include "comsol_text_exporter.hpp"
#include "logger.hpp"

#include <fstream>
#include <iomanip>

namespace mpfem {

void ComsolTextExporter::write(const std::string &filePath, const CoupledFieldResult &result)
{
    Check(!result.coordinates.empty(), "ComsolTextExporter received empty result");
    Check(result.electricPotential.size() == result.coordinates.size(), 
          "ComsolTextExporter field size mismatch");
    Check(result.temperature.size() == result.coordinates.size(), 
          "ComsolTextExporter field size mismatch");
    Check(result.displacement.size() == result.coordinates.size(), 
          "ComsolTextExporter field size mismatch");

    std::ofstream file(filePath);
    Check(file.is_open(), "Cannot open output file: " + filePath);

    file << std::setprecision(16);
    file << "% x y z V T disp\n";
    for (std::size_t i = 0; i < result.coordinates.size(); ++i) {
        file << result.coordinates[i].x << " "
             << result.coordinates[i].y << " "
             << result.coordinates[i].z << " "
             << result.electricPotential[i] << " "
             << result.temperature[i] << " "
             << result.displacement[i] << "\n";
    }
}

} // namespace mpfem