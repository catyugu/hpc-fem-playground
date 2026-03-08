#include "material_xml_reader.hpp"
#include "value_parser.hpp"
#include "logger.hpp"

#include "general/tinyxml2.h"

#include <string>

namespace mpfem {

namespace {

void parsePropertySets(const tinyxml2::XMLElement *materialElement,
                       std::map<std::string, double> &target)
{
    const tinyxml2::XMLElement *propertyGroupElement = materialElement->FirstChildElement("propertyGroup");
    while (propertyGroupElement != nullptr) {
        const tinyxml2::XMLElement *setElement = propertyGroupElement->FirstChildElement("set");
        while (setElement != nullptr) {
            const char *propertyName = setElement->Attribute("name");
            const char *valueText = setElement->Attribute("value");
            if (propertyName != nullptr && valueText != nullptr) {
                double parsedValue = 0.0;
                if (ValueParser::parseFirstNumber(valueText, parsedValue)) {
                    target[propertyName] = parsedValue;
                }
            }
            setElement = setElement->NextSiblingElement("set");
        }
        propertyGroupElement = propertyGroupElement->NextSiblingElement("propertyGroup");
    }

    const tinyxml2::XMLElement *setElement = materialElement->FirstChildElement("set");
    while (setElement != nullptr) {
        const char *propertyName = setElement->Attribute("name");
        const char *valueText = setElement->Attribute("value");
        if (propertyName != nullptr && valueText != nullptr) {
            double parsedValue = 0.0;
            if (ValueParser::parseFirstNumber(valueText, parsedValue)) {
                target[propertyName] = parsedValue;
            }
        }
        setElement = setElement->NextSiblingElement("set");
    }
}

} // namespace

void MaterialXmlReader::readFromFile(const std::string &filePath, PhysicsMaterialDatabase &database)
{
    database.byTag.clear();

    tinyxml2::XMLDocument document;
    const tinyxml2::XMLError loadError = document.LoadFile(filePath.c_str());
    Check(loadError == tinyxml2::XML_SUCCESS, "Failed to parse XML file: " + filePath);

    const tinyxml2::XMLElement *archiveElement = document.FirstChildElement("archive");
    Check(archiveElement != nullptr, "Missing <archive> root node in material XML");
    
    const tinyxml2::XMLElement *modelElement = archiveElement->FirstChildElement("model");
    Check(modelElement != nullptr, "Missing <model> node in material XML");

    const tinyxml2::XMLElement *materialElement = modelElement->FirstChildElement("material");
    while (materialElement != nullptr) {
        MaterialPropertyModel material;
        if (materialElement->Attribute("tag") != nullptr) {
            material.tag = materialElement->Attribute("tag");
        }

        const tinyxml2::XMLElement *labelElement = materialElement->FirstChildElement("label");
        if (labelElement != nullptr && labelElement->Attribute("label") != nullptr) {
            material.label = labelElement->Attribute("label");
        }

        parsePropertySets(materialElement, material.properties);

        // Extract specific properties (names match material.xml)
        auto getProperty = [&](const std::string& name) -> double {
            auto it = material.properties.find(name);
            return it != material.properties.end() ? it->second : 0.0;
        };

        material.rho0 = getProperty("rho0");
        material.alpha = getProperty("alpha");
        material.tref = getProperty("Tref");
        material.electricConductivity = getProperty("electricconductivity");
        material.thermalConductivity = getProperty("thermalconductivity");
        material.youngModulus = getProperty("E");
        material.poissonRatio = getProperty("nu");
        material.thermalExpansion = getProperty("thermalexpansioncoefficient");

        database.byTag[material.tag] = material;
        materialElement = materialElement->NextSiblingElement("material");
    }

    Check(!database.byTag.empty(), "No material blocks found in file: " + filePath);
}

} // namespace mpfem