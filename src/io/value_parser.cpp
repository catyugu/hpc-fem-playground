#include "value_parser.hpp"

#include <cctype>
#include <cstdlib>

namespace mpfem {

bool ValueParser::parseFirstNumber(const std::string &text, double &value)
{
    value = 0.0;

    std::string token;
    for (std::size_t i = 0; i < text.size(); ++i) {
        const char current = text[i];
        const bool isNumericChar = std::isdigit(static_cast<unsigned char>(current)) != 0
                                   || current == '+'
                                   || current == '-'
                                   || current == '.'
                                   || current == 'e'
                                   || current == 'E';
        if (isNumericChar) {
            token.push_back(current);
            continue;
        }

        if (!token.empty()) {
            break;
        }
    }

    if (token.empty()) {
        return false;
    }

    char *endPtr = nullptr;
    value = std::strtod(token.c_str(), &endPtr);
    if (endPtr == token.c_str()) {
        return false;
    }

    return true;
}

} // namespace mpfem
