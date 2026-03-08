#include "material_xml_reader.hpp"

#include <string>

int main(int argc, char **argv)
{
    const int EXPECTED_ARGUMENT_COUNT = 2;
    if (argc < EXPECTED_ARGUMENT_COUNT) {
        return 1;
    }

    mpfem::MaterialDatabase database;
    std::string errorMessage;
    if (!mpfem::MaterialXmlReader::readFromFile(argv[1], database, errorMessage)) {
        return 1;
    }

    if (database.materials.size() != 2) {
        return 1;
    }

    const mpfem::MaterialDefinition &mat1 = database.materials[0];
    if (mat1.tag != "mat1") {
        return 1;
    }

    if (mat1.siProperties.find("rho0") == mat1.siProperties.end()) {
        return 1;
    }
    if (mat1.siProperties.find("alpha") == mat1.siProperties.end()) {
        return 1;
    }
    if (mat1.siProperties.find("Tref") == mat1.siProperties.end()) {
        return 1;
    }

    return 0;
}
