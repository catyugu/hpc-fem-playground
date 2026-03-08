#include "material_xml_reader.hpp"

#include "value_parser.hpp"

#include "general/tinyxml2.h"

#include <string>

namespace mpfem {

namespace {
void parsePropertySets(const tinyxml2::XMLElement *materialElement,
                       std::map<std::string, double> &target,
                       std::string &errorMessage)
{
    const tinyxml2::XMLElement *propertyGroupElement = materialElement->FirstChildElement("propertyGroup");
    while (propertyGroupElement != NULL) {
        const tinyxml2::XMLElement *setElement = propertyGroupElement->FirstChildElement("set");
        while (setElement != NULL) {
            const char *propertyName = setElement->Attribute("name");
            const char *valueText = setElement->Attribute("value");
            if (propertyName != NULL && valueText != NULL) {
                double parsedValue = 0.0;
                std::string parseError;
                if (ValueParser::parseFirstNumber(valueText, parsedValue, parseError)) {
                    target[propertyName] = parsedValue;
                }
            }
            setElement = setElement->NextSiblingElement("set");
        }

        propertyGroupElement = propertyGroupElement->NextSiblingElement("propertyGroup");
    }

    const tinyxml2::XMLElement *setElement = materialElement->FirstChildElement("set");
    while (setElement != NULL) {
        const char *propertyName = setElement->Attribute("name");
        const char *valueText = setElement->Attribute("value");
        if (propertyName != NULL && valueText != NULL) {
            double parsedValue = 0.0;
            std::string parseError;
            if (ValueParser::parseFirstNumber(valueText, parsedValue, parseError)) {
                target[propertyName] = parsedValue;
            }
        }
        setElement = setElement->NextSiblingElement("set");
    }

    errorMessage.clear();
}

} // namespace

bool MaterialXmlReader::readFromFile(const std::string &filePath,
                                     MaterialDatabase &database,
                                     std::string &errorMessage)
{
    database.materials.clear();
    errorMessage.clear();

    tinyxml2::XMLDocument document;
    const tinyxml2::XMLError loadError = document.LoadFile(filePath.c_str());
    if (loadError != tinyxml2::XML_SUCCESS) {
        errorMessage = "Failed to parse XML file: " + filePath;
        return false;
    }

    const tinyxml2::XMLElement *archiveElement = document.FirstChildElement("archive");
    if (archiveElement == NULL) {
        errorMessage = "Missing <archive> root node in material XML";
        return false;
    }
    const tinyxml2::XMLElement *modelElement = archiveElement->FirstChildElement("model");
    if (modelElement == NULL) {
        errorMessage = "Missing <model> node in material XML";
        return false;
    }

    const tinyxml2::XMLElement *materialElement = modelElement->FirstChildElement("material");
    while (materialElement != NULL) {
        MaterialDefinition material;
        if (materialElement->Attribute("tag") != NULL) {
            material.tag = materialElement->Attribute("tag");
        }

        const tinyxml2::XMLElement *labelElement = materialElement->FirstChildElement("label");
        if (labelElement != NULL && labelElement->Attribute("label") != NULL) {
            material.label = labelElement->Attribute("label");
        }

        std::string parseError;
        parsePropertySets(materialElement, material.siProperties, parseError);
        if (!parseError.empty()) {
            errorMessage = "Failed to parse material properties for tag " + material.tag + ": " + parseError;
            return false;
        }

        database.materials.push_back(material);
        materialElement = materialElement->NextSiblingElement("material");
    }

    if (database.materials.empty()) {
        errorMessage = "No material blocks found in file: " + filePath;
        return false;
    }

    return true;
}

} // namespace mpfem
