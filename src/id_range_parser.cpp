#include "id_range_parser.hpp"
#include "logger.hpp"

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

void IdRangeParser::parseIds(const std::string &text, std::set<int> &ids)
{
    ids.clear();

    std::stringstream tokenStream(text);
    std::string token;

    while (std::getline(tokenStream, token, TOKEN_DELIMITER)) {
        const std::string trimmedToken = trim(token);
        Check(!trimmedToken.empty(), "Empty token found in range expression");

        const std::size_t rangePos = trimmedToken.find(RANGE_DELIMITER);
        if (rangePos == std::string::npos) {
            int value = 0;
            Check(parsePositiveInt(trimmedToken, value), 
                  ("Invalid identifier token: " + trimmedToken).c_str());
            ids.insert(value);
            continue;
        }

        const std::string beginToken = trim(trimmedToken.substr(0, rangePos));
        const std::string endToken = trim(trimmedToken.substr(rangePos + 1));
        int beginValue = 0;
        int endValue = 0;
        Check(parsePositiveInt(beginToken, beginValue) && parsePositiveInt(endToken, endValue),
              ("Invalid range token: " + trimmedToken).c_str());
        Check(beginValue <= endValue, 
              ("Descending range is not allowed: " + trimmedToken).c_str());

        for (int value = beginValue; value <= endValue; ++value) {
            ids.insert(value);
        }
    }

    Check(!ids.empty(), "No ids found in expression");
}

} // namespace mpfem