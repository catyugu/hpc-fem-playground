#ifndef MPFEM_ID_RANGE_PARSER_HPP
#define MPFEM_ID_RANGE_PARSER_HPP

#include <set>
#include <string>

namespace mpfem {

/**
 * @brief Parses integer identifiers from comma-separated tokens and ranges.
 */
class IdRangeParser {
public:
    /**
     * @brief Parses expressions like "1-3,8,10-12" into a unique id set.
     * @param text Input expression.
     * @param ids Parsed ids.
     * @param errorMessage Error description when parsing fails.
     * @return True when parsing succeeds.
     */
    static bool parseIds(const std::string &text,
                         std::set<int> &ids,
                         std::string &errorMessage);
};

} // namespace mpfem

#endif
