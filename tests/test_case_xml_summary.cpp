#include "case_xml_summary_reader.hpp"

#include <string>

int main(int argc, char **argv)
{
    const int EXPECTED_ARGUMENT_COUNT = 2;
    if (argc < EXPECTED_ARGUMENT_COUNT) {
        return 1;
    }

    mpfem::CaseXmlSummary summary;
    std::string error;
    if (!mpfem::CaseXmlSummaryReader::readFromFile(argv[1], summary, error)) {
        return 1;
    }

    const int EXPECTED_VARIABLES = 7;
    const std::size_t EXPECTED_DOMAINS = 7;
    const std::size_t EXPECTED_BOUNDARIES = 43;

    if (summary.variableCount != EXPECTED_VARIABLES) {
        return 1;
    }
    if (summary.domainIds.size() != EXPECTED_DOMAINS) {
        return 1;
    }
    if (summary.boundaryIds.size() != EXPECTED_BOUNDARIES) {
        return 1;
    }

    if (summary.domainIds.find(1) == summary.domainIds.end() ||
        summary.domainIds.find(7) == summary.domainIds.end()) {
        return 1;
    }
    if (summary.boundaryIds.find(1) == summary.boundaryIds.end() ||
        summary.boundaryIds.find(43) == summary.boundaryIds.end()) {
        return 1;
    }

    return 0;
}
