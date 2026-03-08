#include "case_xml_reader.hpp"

#include <string>

int main(int argc, char **argv)
{
    const int EXPECTED_ARGUMENT_COUNT = 2;
    if (argc < EXPECTED_ARGUMENT_COUNT) {
        return 1;
    }

    mpfem::CaseDefinition caseDefinition;
    std::string errorMessage;
    if (!mpfem::CaseXmlReader::readFromFile(argv[1], caseDefinition, errorMessage)) {
        return 1;
    }

    if (caseDefinition.caseName != "busbar") {
        return 1;
    }
    if (caseDefinition.variables.size() != 7) {
        return 1;
    }
    if (caseDefinition.materialAssignments.size() != 2) {
        return 1;
    }
    if (caseDefinition.physicsDefinitions.size() != 3) {
        return 1;
    }
    if (caseDefinition.coupledPhysicsDefinitions.size() != 2) {
        return 1;
    }

    return 0;
}
