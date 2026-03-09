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
    mfem::Mpi::Init(argc, argv);
    mfem::Hypre::Init();
#endif

    int rank = 0;
#ifdef MFEM_USE_MPI
    rank = mfem::Mpi::WorldRank();
#endif

    if (argc < 2) {
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
    FemMesh mesh;
    {
        ScopedTimer timer("Mesh loading");

        std::filesystem::path mfemMeshPath = meshPath;
        if (meshPath.extension() == ".mphtxt") {
            mfemMeshPath.replace_extension(".mesh");

#ifdef MFEM_USE_MPI
            if (rank == 0) {
                const std::string cmd = "python3 scripts/mphtxt_to_mfem_mesh.py "
                                        + meshPath.string() + " " + mfemMeshPath.string();
                if (std::system(cmd.c_str()) != 0) {
                    Logger::log(LogLevel::Error, "Failed to convert mphtxt mesh");
                    return 1;
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
#else
            const std::string cmd = "python3 scripts/mphtxt_to_mfem_mesh.py "
                                   + meshPath.string() + " " + mfemMeshPath.string();
            if (std::system(cmd.c_str()) != 0) {
                Logger::log(LogLevel::Error, "Failed to convert mphtxt mesh");
                return 1;
            }
#endif
        }

        auto serialMesh = std::make_unique<mfem::Mesh>(mfemMeshPath.string().c_str(), 1, 1, true);
        
#ifdef MFEM_USE_MPI
        mesh = FemMesh(MPI_COMM_WORLD, *serialMesh);
#else
        mesh = std::move(*serialMesh);
#endif
    }

    // ========== Step 6: Run coupled solver ==========
    MfemCoupledSolver solver;
    CoupledFieldResult result;
    solver.solve(mesh, problemModel, materialDatabase, result);

    // ========== Step 7: Export results ==========
    std::filesystem::path comsolOutput;
    {
        ScopedTimer timer("Result export");

#ifdef MFEM_USE_MPI
        if (rank == 0) {
#endif
            std::filesystem::path resultsDir = std::filesystem::current_path() / "results" / caseDefinition.caseName;
            std::filesystem::create_directories(resultsDir);

            comsolOutput = resultsDir / "mpfem_result.txt";
            ComsolTextExporter::write(comsolOutput.string(), result);

            std::filesystem::path vtuOutput = resultsDir / "mpfem_result.vtu";
            VtuExporter::write(vtuOutput.string(), result);
#ifdef MFEM_USE_MPI
        }
        MPI_Barrier(MPI_COMM_WORLD);
#endif
    }

    if (rank == 0) {
        Logger::log(LogLevel::Info, "Output: " + comsolOutput.string());
    }

    return 0;
}