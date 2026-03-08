#include "case_xml_summary_reader.hpp"
#include "logger.hpp"

#include <string>

int main(int argc, char **argv)
{
    const int EXPECTED_ARGUMENT_COUNT = 2;
    if (argc < EXPECTED_ARGUMENT_COUNT) {
        mpfem::Logger::log(mpfem::LogLevel::Error,
                            "Usage: busbar_case_summary <path-to-case.xml>");
        return 1;
    }

    mpfem::CaseXmlSummary summary;
    std::string errorMessage;
    if (!mpfem::CaseXmlSummaryReader::readFromFile(argv[1], summary, errorMessage)) {
        mpfem::Logger::log(mpfem::LogLevel::Error, errorMessage);
        return 1;
    }

    mpfem::Logger::log(mpfem::LogLevel::Info,
                        "Parsed case summary: variables=" + std::to_string(summary.variableCount)
                            + ", domains=" + std::to_string(summary.domainIds.size())
                            + ", boundaries=" + std::to_string(summary.boundaryIds.size()));
    return 0;
}
