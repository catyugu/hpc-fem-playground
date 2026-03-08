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
     * @param errorMessage Error details on failure.
     * @return True when parsing succeeds.
     */
    static bool parseFirstNumber(const std::string &text,
                                 double &value,
                                 std::string &errorMessage);
};

} // namespace mpfem

#endif
