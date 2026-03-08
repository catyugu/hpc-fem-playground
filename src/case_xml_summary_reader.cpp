#include "case_xml_summary_reader.hpp"

#include "id_range_parser.hpp"

#include <fstream>
#include <regex>
#include <sstream>

namespace mpfem {

namespace {

const char *const VAR_PATTERN = "<var\\s+";
const char *const ASSIGN_DOMAINS_PATTERN = "<assign[^>]*domains=\\\"([^\\\"]+)\\\"";
const char *const BOUNDARY_IDS_PATTERN = "<boundary[^>]*ids=\\\"([^\\\"]+)\\\"";

bool readTextFile(const std::string &filePath,
                  std::string &content,
                  std::string &errorMessage)
{
    std::ifstream input(filePath.c_str());
    if (!input.is_open()) {
        errorMessage = "Cannot open file: " + filePath;
        return false;
    }

    std::ostringstream buffer;
    buffer << input.rdbuf();
    content = buffer.str();
    if (content.empty()) {
        errorMessage = "Empty file: " + filePath;
        return false;
    }

    return true;
}

bool collectIdsByPattern(const std::string &content,
                         const std::string &pattern,
                         std::set<int> &target,
                         std::string &errorMessage)
{
    const std::regex regexPattern(pattern);
    const std::sregex_iterator begin(content.begin(), content.end(), regexPattern);
    const std::sregex_iterator end;

    for (std::sregex_iterator iter = begin; iter != end; ++iter) {
        const std::string rangeText = (*iter)[1].str();
        std::set<int> parsedIds;
        std::string parseError;
        if (!IdRangeParser::parseIds(rangeText, parsedIds, parseError)) {
            errorMessage = "Failed to parse id list '" + rangeText + "': " + parseError;
            return false;
        }

        target.insert(parsedIds.begin(), parsedIds.end());
    }

    return true;
}

} // namespace

bool CaseXmlSummaryReader::readFromFile(const std::string &filePath,
                                        CaseXmlSummary &summary,
                                        std::string &errorMessage)
{
    summary.variableCount = 0;
    summary.domainIds.clear();
    summary.boundaryIds.clear();
    errorMessage.clear();

    std::string content;
    if (!readTextFile(filePath, content, errorMessage)) {
        return false;
    }

    const std::regex variablePattern(VAR_PATTERN);
    const std::sregex_iterator begin(content.begin(), content.end(), variablePattern);
    const std::sregex_iterator end;
    int count = 0;
    for (std::sregex_iterator iter = begin; iter != end; ++iter) {
        (void)iter;
        ++count;
    }
    summary.variableCount = count;

    if (!collectIdsByPattern(content,
                             ASSIGN_DOMAINS_PATTERN,
                             summary.domainIds,
                             errorMessage)) {
        return false;
    }
    if (!collectIdsByPattern(content,
                             BOUNDARY_IDS_PATTERN,
                             summary.boundaryIds,
                             errorMessage)) {
        return false;
    }

    return true;
}

} // namespace mpfem
