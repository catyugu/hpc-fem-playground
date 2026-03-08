#include "case_xml_reader.hpp"
#include "comsol_result_reader.hpp"
#include "comsol_text_exporter.hpp"
#include "logger.hpp"
#include "material_xml_reader.hpp"
#include "mfem_coupled_solver.hpp"
#include "physics_problem_builder.hpp"
#include "vtu_exporter.hpp"

#include "mfem.hpp"

#include <filesystem>
#include <cstdlib>
#include <string>
#include <vector>

int main(int argc, char **argv)
{
    const int EXPECTED_ARGUMENT_COUNT = 2;
    if (argc < EXPECTED_ARGUMENT_COUNT) {
        mpfem::Logger::log(mpfem::LogLevel::Error,
                           "Usage: busbar_pipeline <path-to-case.xml>");
        return 1;
    }

    const std::filesystem::path casePath(argv[1]);
    const std::filesystem::path caseDirectory = casePath.parent_path();

    mpfem::CaseDefinition caseDefinition;
    std::string errorMessage;
    if (!mpfem::CaseXmlReader::readFromFile(casePath.string(), caseDefinition, errorMessage)) {
        mpfem::Logger::log(mpfem::LogLevel::Error, errorMessage);
        return 1;
    }

    const std::filesystem::path materialPath = caseDirectory / caseDefinition.materialsPath;
    const std::filesystem::path referencePath = caseDirectory / caseDefinition.comsolResultPath;
    const std::filesystem::path meshPath = caseDirectory / caseDefinition.meshPath;

    mpfem::MaterialDatabase materialDatabase;
    if (!mpfem::MaterialXmlReader::readFromFile(materialPath.string(), materialDatabase, errorMessage)) {
        mpfem::Logger::log(mpfem::LogLevel::Error, errorMessage);
        return 1;
    }

    mpfem::PhysicsProblemModel problemModel;
    mpfem::PhysicsMaterialDatabase physicsMaterials;
    if (!mpfem::PhysicsProblemBuilder::build(caseDefinition,
                                             materialDatabase,
                                             problemModel,
                                             physicsMaterials,
                                             errorMessage)) {
        mpfem::Logger::log(mpfem::LogLevel::Error, errorMessage);
        return 1;
    }

    std::vector<mpfem::ComsolResultRow> referenceRows;
    if (!mpfem::ComsolResultReader::readFromFile(referencePath.string(), referenceRows, errorMessage)) {
        mpfem::Logger::log(mpfem::LogLevel::Error, errorMessage);
        return 1;
    }

    std::filesystem::path mfemMeshPath = meshPath;
    if (meshPath.extension() == ".mphtxt") {
        mfemMeshPath = meshPath;
        mfemMeshPath.replace_extension(".msh");
        const std::string convertCommand = std::string("python3 scripts/mphtxt_to_gmsh.py ")
                                           + meshPath.string()
                                           + " "
                                           + mfemMeshPath.string();
        const int convertStatus = std::system(convertCommand.c_str());
        if (convertStatus != 0) {
            mpfem::Logger::log(mpfem::LogLevel::Error,
                               "Failed to convert mphtxt mesh to gmsh format");
            return 1;
        }
    }

    mfem::Mesh mesh(mfemMeshPath.string().c_str(), 1, 1, true);

    mpfem::MfemCoupledSolver solver;
    mpfem::CoupledFieldResult result;
    if (!solver.solve(mesh,
                      problemModel,
                      physicsMaterials,
                      12,
                      1e-7,
                      result,
                      errorMessage)) {
        mpfem::Logger::log(mpfem::LogLevel::Error, errorMessage);
        return 1;
    }

    const std::filesystem::path outputDirectory = std::filesystem::path("results") / caseDefinition.caseName;
    std::filesystem::create_directories(outputDirectory);

    const std::filesystem::path comsolOutput = outputDirectory / "mpfem_result.txt";
    if (!mpfem::ComsolTextExporter::write(comsolOutput.string(), result, errorMessage)) {
        mpfem::Logger::log(mpfem::LogLevel::Error, errorMessage);
        return 1;
    }

    const std::filesystem::path vtuOutput = outputDirectory / "mpfem_result.vtu";
    if (!mpfem::VtuExporter::write(vtuOutput.string(), result, errorMessage)) {
        mpfem::Logger::log(mpfem::LogLevel::Error, errorMessage);
        return 1;
    }

    mpfem::Logger::log(mpfem::LogLevel::Info,
                       "Busbar baseline pipeline completed. Output: " + comsolOutput.string());
    return 0;
}
