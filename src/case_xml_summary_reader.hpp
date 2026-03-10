#ifndef MPFEM_CASE_XML_SUMMARY_READER_HPP
#define MPFEM_CASE_XML_SUMMARY_READER_HPP

#include <set>
#include <string>

namespace mpfem {

/**
 * @brief Minimal extracted metadata from case XML used by early tests.
 */
struct CaseXmlSummary {
    int variableCount = 0;
    std::set<int> domainIds;
    std::set<int> boundaryIds;
};

/**
 * @brief Reads lightweight counts and id sets from case XML.
 */
class CaseXmlSummaryReader {
public:
    /**
     * @brief Loads summary data from XML text.
     */
    static void readFromFile(const std::string &filePath, CaseXmlSummary &summary);
};

} // namespace mpfem

#endif