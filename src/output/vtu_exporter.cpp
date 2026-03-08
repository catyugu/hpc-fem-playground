#include "vtu_exporter.hpp"
#include "logger.hpp"

#include <fstream>

namespace mpfem {

void VtuExporter::write(const std::string &filePath, const CoupledFieldResult &result)
{
    Check(!result.coordinates.empty(), "VtuExporter received empty result");
    Check(result.electricPotential.size() == result.coordinates.size(), 
          "VtuExporter field size mismatch");
    Check(result.temperature.size() == result.coordinates.size(), 
          "VtuExporter field size mismatch");
    Check(result.displacement.size() == result.coordinates.size(), 
          "VtuExporter field size mismatch");

    std::ofstream file(filePath);
    Check(file.is_open(), "Cannot open output file: " + filePath);

    const std::size_t numPoints = result.coordinates.size();

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
    file << "<UnstructuredGrid>\n";
    file << "<Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"0\">\n";

    file << "<PointData>\n";
    file << "<DataArray type=\"Float64\" Name=\"ElectricPotential\" format=\"ascii\">\n";
    for (const auto& v : result.electricPotential) {
        file << v << " ";
    }
    file << "\n</DataArray>\n";
    file << "<DataArray type=\"Float64\" Name=\"Temperature\" format=\"ascii\">\n";
    for (const auto& t : result.temperature) {
        file << t << " ";
    }
    file << "\n</DataArray>\n";
    file << "<DataArray type=\"Float64\" Name=\"Displacement\" format=\"ascii\">\n";
    for (const auto& d : result.displacement) {
        file << d << " ";
    }
    file << "\n</DataArray>\n";
    file << "</PointData>\n";

    file << "<Points>\n";
    file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto& coord : result.coordinates) {
        file << coord.x << " " << coord.y << " " << coord.z << " ";
    }
    file << "\n</DataArray>\n";
    file << "</Points>\n";

    file << "<Cells>\n";
    file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    file << "</DataArray>\n";
    file << "<DataArray type=\"Int32\" Name=\"types\" format=\"ascii\">\n";
    file << "</DataArray>\n";
    file << "</Cells>\n";

    file << "</Piece>\n";
    file << "</UnstructuredGrid>\n";
    file << "</VTKFile>\n";
}

} // namespace mpfem