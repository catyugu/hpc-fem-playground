#include "vtu_exporter.hpp"

#include <fstream>

namespace mpfem {

bool VtuExporter::write(const std::string &outputPath,
                        const CoupledFieldResult &result,
                        std::string &errorMessage)
{
    errorMessage.clear();

    const std::size_t count = result.coordinates.size();
    if (count == 0) {
        errorMessage = "VtuExporter received empty result";
        return false;
    }
    if (result.electricPotential.size() != count
        || result.temperature.size() != count
        || result.displacement.size() != count) {
        errorMessage = "VtuExporter field size mismatch";
        return false;
    }

    std::ofstream output(outputPath.c_str());
    if (!output.is_open()) {
        errorMessage = "Cannot open output file: " + outputPath;
        return false;
    }

    output << "<?xml version=\"1.0\"?>\n";
    output << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    output << "  <UnstructuredGrid>\n";
    output << "    <Piece NumberOfPoints=\"" << count << "\" NumberOfCells=\"" << count << "\">\n";

    output << "      <Points>\n";
    output << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < count; ++i) {
        output << "          " << result.coordinates[i].x << " "
               << result.coordinates[i].y << " "
               << result.coordinates[i].z << "\n";
    }
    output << "        </DataArray>\n";
    output << "      </Points>\n";

    output << "      <Cells>\n";
    output << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < count; ++i) {
        output << "          " << i << "\n";
    }
    output << "        </DataArray>\n";

    output << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < count; ++i) {
        output << "          " << (i + 1) << "\n";
    }
    output << "        </DataArray>\n";

    output << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < count; ++i) {
        output << "          1\n";
    }
    output << "        </DataArray>\n";
    output << "      </Cells>\n";

    output << "      <PointData Scalars=\"ElectricPotential Temperature Displacement\">\n";

    output << "        <DataArray type=\"Float64\" Name=\"ElectricPotential\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < count; ++i) {
        output << "          " << result.electricPotential[i] << "\n";
    }
    output << "        </DataArray>\n";

    output << "        <DataArray type=\"Float64\" Name=\"Temperature\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < count; ++i) {
        output << "          " << result.temperature[i] << "\n";
    }
    output << "        </DataArray>\n";

    output << "        <DataArray type=\"Float64\" Name=\"Displacement\" format=\"ascii\">\n";
    for (std::size_t i = 0; i < count; ++i) {
        output << "          " << result.displacement[i] << "\n";
    }
    output << "        </DataArray>\n";

    output << "      </PointData>\n";
    output << "    </Piece>\n";
    output << "  </UnstructuredGrid>\n";
    output << "</VTKFile>\n";

    return true;
}

} // namespace mpfem
