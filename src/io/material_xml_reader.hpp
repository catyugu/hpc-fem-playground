#ifndef MPFEM_MATERIAL_XML_READER_HPP
#define MPFEM_MATERIAL_XML_READER_HPP

#include "material_model.hpp"

#include <string>

namespace mpfem {

/**
 * @brief Reads material XML and extracts SI scalar values.
 */
class MaterialXmlReader {
public:
    /**
     * @brief Parses material XML from disk.
     * @param filePath Input material XML path.
     * @param database Parsed material database.
     * @param errorMessage Error details on failure.
     * @return True when parsing succeeds.
     */
    static bool readFromFile(const std::string &filePath,
                             MaterialDatabase &database,
                             std::string &errorMessage);
};

} // namespace mpfem

#endif
