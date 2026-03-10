#ifndef MPFEM_MATERIAL_XML_READER_HPP
#define MPFEM_MATERIAL_XML_READER_HPP

#include "model/physics_problem_model.hpp"

#include <string>

namespace mpfem {

/**
 * @brief Reads material properties from material XML.
 */
class MaterialXmlReader {
public:
    /**
     * @brief Parses material XML from disk.
     */
    static void readFromFile(const std::string &filePath, PhysicsMaterialDatabase &database);
};

} // namespace mpfem

#endif