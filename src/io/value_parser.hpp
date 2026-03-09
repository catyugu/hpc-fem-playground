#ifndef MPFEM_VALUE_PARSER_HPP
#define MPFEM_VALUE_PARSER_HPP

#include <string>

namespace mpfem {

/**
 * @brief Parses numeric prefix from COMSOL-like value strings.
 */
class ValueParser {
public:
    /**
     * @brief Parses first scalar number from text such as "20[mV]".
     * @param text Input text.
     * @param value Parsed numeric value.
     * @return True when parsing succeeds, false otherwise.
     */
    static bool parseFirstNumber(const std::string &text, double &value);
};

} // namespace mpfem

#endif
