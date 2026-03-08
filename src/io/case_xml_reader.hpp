#ifndef MPFEM_CASE_XML_READER_HPP
#define MPFEM_CASE_XML_READER_HPP

#include "case_model.hpp"

#include <string>

namespace mpfem {

/**
 * @brief Reads complete case schema from case XML.
 */
class CaseXmlReader {
public:
    /**
     * @brief Parses case XML from disk.
     * @param filePath Input case XML path.
     * @param caseDefinition Parsed output model.
     * @param errorMessage Error details on failure.
     * @return True when parsing succeeds.
     */
    static bool readFromFile(const std::string &filePath,
                             CaseDefinition &caseDefinition,
                             std::string &errorMessage);
};

} // namespace mpfem

#endif
