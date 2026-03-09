#include "case_xml_summary_reader.hpp"
#include "id_range_parser.hpp"
#include "logger.hpp"

#include <fstream>
#include <regex>
#include <sstream>

namespace mpfem {

namespace {

const char *const VAR_PATTERN = "<var\\s+";
const char *const ASSIGN_DOMAINS_PATTERN = "<assign[^>]*domains=\\\"([^\\\"]+)\\\"";
const char *const BOUNDARY_IDS_PATTERN = "<boundary[^>]*ids=\\\"([^\\\"]+)\\\"";

void readTextFile(const std::string &filePath, std::string &content)
{
    std::ifstream input(filePath.c_str());
    Check(input.is_open(), "Cannot open file: " + filePath);

    std::ostringstream buffer;
    buffer << input.rdbuf();
    content = buffer.str();
    Check(!content.empty(), "Empty file: " + filePath);
}

void collectIdsByPattern(const std::string &content,
                         const std::string &pattern,
                         std::set<int> &target)
{
    const std::regex regexPattern(pattern);
    const std::sregex_iterator begin(content.begin(), content.end(), regexPattern);
    const std::sregex_iterator end;

    for (std::sregex_iterator iter = begin; iter != end; ++iter) {
        const std::string rangeText = (*iter)[1].str();
        std::set<int> parsedIds;
        IdRangeParser::parseIds(rangeText, parsedIds);
        target.insert(parsedIds.begin(), parsedIds.end());
    }
}

} // namespace

void CaseXmlSummaryReader::readFromFile(const std::string &filePath, CaseXmlSummary &summary)
{
    summary.variableCount = 0;
    summary.domainIds.clear();
    summary.boundaryIds.clear();

    std::string content;
    readTextFile(filePath, content);

    const std::regex variablePattern(VAR_PATTERN);
    const std::sregex_iterator begin(content.begin(), content.end(), variablePattern);
    const std::sregex_iterator end;
    int count = 0;
    for (std::sregex_iterator iter = begin; iter != end; ++iter) {
        (void)iter;
        ++count;
    }
    summary.variableCount = count;

    collectIdsByPattern(content, ASSIGN_DOMAINS_PATTERN, summary.domainIds);
    collectIdsByPattern(content, BOUNDARY_IDS_PATTERN, summary.boundaryIds);
}

} // namespace mpfem