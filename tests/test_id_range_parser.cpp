#include "id_range_parser.hpp"

#include <set>
#include <string>

int main()
{
    std::set<int> ids;
    std::string error;

    if (!mpfem::IdRangeParser::parseIds("1-3,8,10-12", ids, error)) {
        return 1;
    }

    const std::size_t EXPECTED_COUNT = 7;
    if (ids.size() != EXPECTED_COUNT) {
        return 1;
    }

    if (ids.find(1) == ids.end() || ids.find(12) == ids.end()) {
        return 1;
    }

    if (mpfem::IdRangeParser::parseIds("3-1", ids, error)) {
        return 1;
    }

    return 0;
}
