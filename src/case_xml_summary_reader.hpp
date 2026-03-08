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
     * @param filePath Input XML file path.
     * @param summary Output summary.
     * @param errorMessage Error details when read fails.
     * @return True when summary is available.
     */
    static bool readFromFile(const std::string &filePath,
                             CaseXmlSummary &summary,
                             std::string &errorMessage);
};

} // namespace mpfem

#endif
