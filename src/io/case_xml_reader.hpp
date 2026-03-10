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
     */
    static void readFromFile(const std::string &filePath, CaseDefinition &caseDefinition);
};

} // namespace mpfem

#endif