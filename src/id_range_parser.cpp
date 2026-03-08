#include "id_range_parser.hpp"

#include <algorithm>
#include <cctype>
#include <sstream>

namespace mpfem {

namespace {

const char TOKEN_DELIMITER = ',';
const char RANGE_DELIMITER = '-';
const int MIN_VALID_ID = 1;

std::string trim(const std::string &text)
{
    std::size_t first = 0;
    while (first < text.size() && std::isspace(static_cast<unsigned char>(text[first])) != 0) {
        ++first;
    }

    std::size_t last = text.size();
    while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1])) != 0) {
        --last;
    }

    return text.substr(first, last - first);
}

bool parsePositiveInt(const std::string &text, int &value)
{
    if (text.empty()) {
        return false;
    }

    for (std::size_t i = 0; i < text.size(); ++i) {
        if (std::isdigit(static_cast<unsigned char>(text[i])) == 0) {
            return false;
        }
    }

    try {
        value = std::stoi(text);
    } catch (...) {
        return false;
    }

    return value >= MIN_VALID_ID;
}

} // namespace

bool IdRangeParser::parseIds(const std::string &text,
                             std::set<int> &ids,
                             std::string &errorMessage)
{
    ids.clear();
    errorMessage.clear();

    std::stringstream tokenStream(text);
    std::string token;

    while (std::getline(tokenStream, token, TOKEN_DELIMITER)) {
        const std::string trimmedToken = trim(token);
        if (trimmedToken.empty()) {
            errorMessage = "Empty token found in range expression";
            return false;
        }

        const std::size_t rangePos = trimmedToken.find(RANGE_DELIMITER);
        if (rangePos == std::string::npos) {
            int value = 0;
            if (!parsePositiveInt(trimmedToken, value)) {
                errorMessage = "Invalid identifier token: " + trimmedToken;
                return false;
            }
            ids.insert(value);
            continue;
        }

        const std::string beginToken = trim(trimmedToken.substr(0, rangePos));
        const std::string endToken = trim(trimmedToken.substr(rangePos + 1));
        int beginValue = 0;
        int endValue = 0;
        if (!parsePositiveInt(beginToken, beginValue) || !parsePositiveInt(endToken, endValue)) {
            errorMessage = "Invalid range token: " + trimmedToken;
            return false;
        }
        if (beginValue > endValue) {
            errorMessage = "Descending range is not allowed: " + trimmedToken;
            return false;
        }

        for (int value = beginValue; value <= endValue; ++value) {
            ids.insert(value);
        }
    }

    if (ids.empty()) {
        errorMessage = "No ids found in expression";
        return false;
    }

    return true;
}

} // namespace mpfem
