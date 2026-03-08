#include "comsol_text_exporter.hpp"

#include <fstream>
#include <iomanip>

namespace mpfem {

bool ComsolTextExporter::write(const std::string &outputPath,
                               const CoupledFieldResult &result,
                               std::string &errorMessage)
{
    errorMessage.clear();

    const std::size_t count = result.coordinates.size();
    if (count == 0) {
        errorMessage = "ComsolTextExporter received empty result";
        return false;
    }
    if (result.electricPotential.size() != count
        || result.temperature.size() != count
        || result.displacement.size() != count) {
        errorMessage = "ComsolTextExporter field size mismatch";
        return false;
    }

    std::ofstream output(outputPath.c_str());
    if (!output.is_open()) {
        errorMessage = "Cannot open output file: " + outputPath;
        return false;
    }

    output << "% Model: mpfem\n";
    output << "% Dimension: 3\n";
    output << "% Nodes: " << count << "\n";
    output << "% Expressions: 3\n";
    output << "% Description: Electric potential, Temperature, Displacement magnitude\n";
    output << "x\ty\tz\tV(V)\tT(K)\tsolid.disp(m)\n";

    output << std::setprecision(16);
    for (std::size_t i = 0; i < count; ++i) {
        output << result.coordinates[i].x << '\t'
               << result.coordinates[i].y << '\t'
               << result.coordinates[i].z << '\t'
               << result.electricPotential[i] << '\t'
               << result.temperature[i] << '\t'
               << result.displacement[i] << '\n';
    }

    return true;
}

} // namespace mpfem
