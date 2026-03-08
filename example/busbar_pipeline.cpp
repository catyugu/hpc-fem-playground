/**
 * @file busbar_pipeline.cpp
 * @brief Complete pipeline for the busbar electro-thermal-mechanical analysis.
 */

#include "mpfem.hpp"
#include "mpfem_types.hpp"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <cstdlib>

using namespace mpfem;

int main(int argc, char* argv[])
{
#ifdef MFEM_USE_MPI
    // Initialize MPI and HYPRE
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();
#endif

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <case.xml>" << std::endl;
        return 1;
    }

    const std::filesystem::path casePath(argv[1]);
    const std::filesystem::path caseDir = casePath.parent_path();

    // ========== Step 1: Parse case XML ==========
    ScopedTimer xmlTimer("XML parsing");
    CaseDefinition caseDefinition;
    CaseXmlReader::readFromFile(casePath.string(), caseDefinition);
    xmlTimer.stop();
    Logger::log(LogLevel::Info, "XML parsing completed in " + xmlTimer.getElapsedStr());

    // Use paths from case XML
    std::filesystem::path materialPath = caseDir / caseDefinition.materialsPath;
    std::filesystem::path meshPath = caseDir / caseDefinition.meshPath;
    std::filesystem::path referencePath = caseDir / caseDefinition.comsolResultPath;

    // ========== Step 2: Load materials ==========
    ScopedTimer materialTimer("Material loading");
    PhysicsMaterialDatabase materialDatabase;
    MaterialXmlReader::readFromFile(materialPath.string(), materialDatabase);
    materialTimer.stop();
    Logger::log(LogLevel::Info, "Material loading completed in " + materialTimer.getElapsedStr());

    // ========== Step 3: Build problem model ==========
    ScopedTimer buildTimer("Problem building");
    PhysicsProblemModel problemModel;
    PhysicsProblemBuilder::build(caseDefinition, materialDatabase, problemModel);
    buildTimer.stop();
    Logger::log(LogLevel::Info, "Problem building completed in " + buildTimer.getElapsedStr());

    // ========== Step 4: Load reference result ==========
    ScopedTimer refTimer("Reference result loading");
    std::vector<ComsolResultRow> referenceRows;
    ComsolResultReader::readFromFile(referencePath.string(), referenceRows);
    refTimer.stop();
    Logger::log(LogLevel::Info, "Reference result loading completed in " + refTimer.getElapsedStr());

    // ========== Step 5: Load mesh ==========
    ScopedTimer meshTimer("Mesh loading");
    
    // Convert mphtxt to MFEM mesh format if needed
    std::filesystem::path mfemMeshPath = meshPath;
    if (meshPath.extension() == ".mphtxt") {
        mfemMeshPath = meshPath;
        mfemMeshPath.replace_extension(".mesh");
        const std::string convertCommand = std::string("python3 scripts/mphtxt_to_mfem_mesh.py ")
                                           + meshPath.string()
                                           + " "
                                           + mfemMeshPath.string();
        const int convertStatus = std::system(convertCommand.c_str());
        if (convertStatus != 0) {
            Logger::log(LogLevel::Error, "Failed to convert mphtxt mesh to MFEM format");
            return 1;
        }
    }

    std::unique_ptr<mfem::Mesh> serialMesh;
    serialMesh.reset(new mfem::Mesh(mfemMeshPath.string().c_str(), 1, 1, true));
    
    FemMesh mesh;
#ifdef MFEM_USE_MPI
    // Create parallel mesh
    mesh = FemMesh(MPI_COMM_WORLD, *serialMesh);
    serialMesh.reset();
#else
    // Use serial mesh directly
    mesh = std::move(*serialMesh);
#endif
    meshTimer.stop();
    Logger::log(LogLevel::Info, "Mesh loading completed in " + meshTimer.getElapsedStr());

    // ========== Step 6: Run coupled solver ==========
    MfemCoupledSolver solver;
    CoupledFieldResult result;
    solver.solve(mesh, problemModel, materialDatabase, result);

    // ========== Step 7: Export results ==========
    const std::filesystem::path resultsDir = std::filesystem::current_path() / "results" / caseDefinition.caseName;
    std::filesystem::create_directories(resultsDir);
    
    ScopedTimer exportTimer("Result export");
    const std::filesystem::path comsolOutput = resultsDir / "mpfem_result.txt";
    ComsolTextExporter::write(comsolOutput.string(), result);
    
    const std::filesystem::path vtuOutput = resultsDir / "mpfem_result.vtu";
    VtuExporter::write(vtuOutput.string(), result);
    
    exportTimer.stop();
    Logger::log(LogLevel::Info, "Result export completed in " + exportTimer.getElapsedStr());
    Logger::log(LogLevel::Info, "Output: " + comsolOutput.string());

    return 0;
}
