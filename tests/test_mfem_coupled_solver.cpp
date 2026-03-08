#include "case_xml_reader.hpp"
#include "material_xml_reader.hpp"
#include "mfem_coupled_solver.hpp"
#include "physics_problem_builder.hpp"

#include "mfem.hpp"

#include <cstdlib>
#include <filesystem>
#include <string>

int main(int argc, char **argv)
{
    const int EXPECTED_ARGUMENT_COUNT = 3;
    if (argc < EXPECTED_ARGUMENT_COUNT) {
        return 1;
    }

    const std::filesystem::path casePath(argv[1]);
    const std::filesystem::path materialPath(argv[2]);
    const std::filesystem::path caseDirectory = casePath.parent_path();

    mpfem::CaseDefinition caseDefinition;
    mpfem::MaterialDatabase materialDatabase;
    std::string errorMessage;

    if (!mpfem::CaseXmlReader::readFromFile(casePath.string(), caseDefinition, errorMessage)) {
        return 1;
    }
    if (!mpfem::MaterialXmlReader::readFromFile(materialPath.string(), materialDatabase, errorMessage)) {
        return 1;
    }

    mpfem::PhysicsProblemModel problemModel;
    mpfem::PhysicsMaterialDatabase physicsMaterials;
    if (!mpfem::PhysicsProblemBuilder::build(caseDefinition,
                                             materialDatabase,
                                             problemModel,
                                             physicsMaterials,
                                             errorMessage)) {
        return 1;
    }

    const std::filesystem::path mphtxtPath = caseDirectory / caseDefinition.meshPath;
    std::filesystem::path gmshPath = mphtxtPath;
    gmshPath.replace_extension(".msh");

    const std::string command = std::string("python3 scripts/mphtxt_to_gmsh.py ")
                                + mphtxtPath.string()
                                + " "
                                + gmshPath.string();
    if (std::system(command.c_str()) != 0) {
        return 1;
    }

    mfem::Mesh mesh(gmshPath.string().c_str(), 1, 1, true);

    mpfem::MfemCoupledSolver solver;
    mpfem::CoupledFieldResult result;
    if (!solver.solve(mesh, problemModel, physicsMaterials, 6, 1e-6, result, errorMessage)) {
        return 1;
    }

    if (result.coordinates.empty() || result.electricPotential.empty()
        || result.temperature.empty() || result.displacement.empty()) {
        return 1;
    }

    return 0;
}
