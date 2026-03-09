/**
 * @file busbar_pipeline.cpp
 * @brief Complete pipeline for the busbar electro-thermal-mechanical analysis.
 */

#include "mpfem.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace mpfem;

int main(int argc, char *argv[])
{
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " <case.xml>" << std::endl;
        return 1;
    }

    const std::filesystem::path casePath(argv[1]);
    const std::filesystem::path caseDir = casePath.parent_path();

    // ========== Step 1: Parse case XML ==========
    CaseDefinition caseDefinition;
    {
        ScopedTimer timer("XML parsing");
        CaseXmlReader::readFromFile(casePath.string(), caseDefinition);
    }

    std::filesystem::path materialPath = caseDir / caseDefinition.materialsPath;
    std::filesystem::path meshPath = caseDir / caseDefinition.meshPath;
    std::filesystem::path referencePath = caseDir / caseDefinition.comsolResultPath;

    // ========== Step 2: Load materials ==========
    PhysicsMaterialDatabase materialDatabase;
    {
        ScopedTimer timer("Material loading");
        MaterialXmlReader::readFromFile(materialPath.string(), materialDatabase);
    }

    // ========== Step 3: Build problem model ==========
    PhysicsProblemModel problemModel;
    {
        ScopedTimer timer("Problem building");
        PhysicsProblemBuilder::build(caseDefinition, materialDatabase, problemModel);
    }

    // ========== Step 4: Load reference result ==========
    {
        ScopedTimer timer("Reference result loading");
        std::vector<ComsolResultRow> referenceRows;
        ComsolResultReader::readFromFile(referencePath.string(), referenceRows);
    }

    // ========== Step 5: Load mesh ==========
    mfem::Mesh mesh;
    {
        ScopedTimer timer("Mesh loading");

        std::filesystem::path MeshPath = meshPath;
        if (meshPath.extension() == ".mphtxt")
        {
            MeshPath.replace_extension(".mesh");

            const std::string cmd = "python3 scripts/mphtxt_to_mfem_mesh.py " + meshPath.string() + " " + MeshPath.string();
            if (std::system(cmd.c_str()) != 0)
            {
                Logger::log(LogLevel::Error, "Failed to convert mphtxt mesh");
                return 1;
            }
        }

        auto mesh = std::make_unique<mfem::Mesh>(MeshPath.string().c_str(), 1, 1, true);
    }

    // ========== Step 6: Run coupled solver ==========
    MfemCoupledSolver solver;
    CoupledFieldResult result;
    solver.solve(mesh, problemModel, materialDatabase, result);

    // ========== Step 7: Export results ==========
    std::filesystem::path comsolOutput;
    {
        ScopedTimer timer("Result export");

        std::filesystem::path resultsDir = std::filesystem::current_path() / "results" / caseDefinition.caseName;
        std::filesystem::create_directories(resultsDir);

        comsolOutput = resultsDir / "mpfem_result.txt";
        ComsolTextExporter::write(comsolOutput.string(), result);

        std::filesystem::path vtuOutput = resultsDir / "mpfem_result.vtu";
        VtuExporter::write(vtuOutput.string(), result);
    }

    return 0;
}