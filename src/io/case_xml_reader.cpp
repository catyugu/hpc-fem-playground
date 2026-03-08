#include "case_xml_reader.hpp"

#include "id_range_parser.hpp"
#include "value_parser.hpp"

#include "general/tinyxml2.h"

#include <sstream>
#include <string>

namespace mpfem {

namespace {

bool parseIdAttribute(const char *attributeText,
                      std::set<int> &target,
                      std::string &errorMessage)
{
    if (attributeText == NULL) {
        target.clear();
        return true;
    }

    return IdRangeParser::parseIds(attributeText, target, errorMessage);
}

void splitCsv(const std::string &csv, std::vector<std::string> &tokens)
{
    tokens.clear();
    std::stringstream stream(csv);
    std::string item;
    while (std::getline(stream, item, ',')) {
        std::string trimmed = item;
        while (!trimmed.empty() && (trimmed.front() == ' ' || trimmed.front() == '\t')) {
            trimmed.erase(trimmed.begin());
        }
        while (!trimmed.empty() && (trimmed.back() == ' ' || trimmed.back() == '\t')) {
            trimmed.pop_back();
        }
        if (!trimmed.empty()) {
            tokens.push_back(trimmed);
        }
    }
}

} // namespace

bool CaseXmlReader::readFromFile(const std::string &filePath,
                                 CaseDefinition &caseDefinition,
                                 std::string &errorMessage)
{
    caseDefinition = CaseDefinition();
    errorMessage.clear();

    tinyxml2::XMLDocument document;
    const tinyxml2::XMLError loadError = document.LoadFile(filePath.c_str());
    if (loadError != tinyxml2::XML_SUCCESS) {
        errorMessage = "Failed to parse XML file: " + filePath;
        return false;
    }

    const tinyxml2::XMLElement *caseElement = document.FirstChildElement("case");
    if (caseElement == NULL) {
        errorMessage = "Missing <case> root node";
        return false;
    }
    if (caseElement->Attribute("name") != NULL) {
        caseDefinition.caseName = caseElement->Attribute("name");
    }

    const tinyxml2::XMLElement *studyElement = caseElement->FirstChildElement("study");
    if (studyElement != NULL && studyElement->Attribute("type") != NULL) {
        caseDefinition.studyType = studyElement->Attribute("type");
    }

    const tinyxml2::XMLElement *pathsElement = caseElement->FirstChildElement("paths");
    if (pathsElement != NULL) {
        if (pathsElement->Attribute("mesh") != NULL) {
            caseDefinition.meshPath = pathsElement->Attribute("mesh");
        }
        if (pathsElement->Attribute("materials") != NULL) {
            caseDefinition.materialsPath = pathsElement->Attribute("materials");
        }
        if (pathsElement->Attribute("comsol_result") != NULL) {
            caseDefinition.comsolResultPath = pathsElement->Attribute("comsol_result");
        }
    }

    const tinyxml2::XMLElement *variablesElement = caseElement->FirstChildElement("variables");
    if (variablesElement != NULL) {
        const tinyxml2::XMLElement *variableElement = variablesElement->FirstChildElement("var");
        while (variableElement != NULL) {
            VariableEntry entry;
            if (variableElement->Attribute("name") != NULL) {
                entry.name = variableElement->Attribute("name");
            }
            if (variableElement->Attribute("value") != NULL) {
                entry.valueText = variableElement->Attribute("value");
            }
            if (variableElement->Attribute("si") != NULL) {
                double parsedSi = 0.0;
                std::string parseError;
                if (!ValueParser::parseFirstNumber(variableElement->Attribute("si"), parsedSi, parseError)) {
                    errorMessage = "Failed to parse variable SI value for " + entry.name + ": " + parseError;
                    return false;
                }
                entry.siValue = parsedSi;
            }
            caseDefinition.variables.push_back(entry);
            variableElement = variableElement->NextSiblingElement("var");
        }
    }

    const tinyxml2::XMLElement *materialsElement = caseElement->FirstChildElement("materials");
    if (materialsElement != NULL) {
        const tinyxml2::XMLElement *assignElement = materialsElement->FirstChildElement("assign");
        while (assignElement != NULL) {
            MaterialAssignment assignment;
            if (assignElement->Attribute("material") != NULL) {
                assignment.materialTag = assignElement->Attribute("material");
            }
            if (!parseIdAttribute(assignElement->Attribute("domains"), assignment.domainIds, errorMessage)) {
                errorMessage = "Failed to parse material assignment domains: " + errorMessage;
                return false;
            }
            caseDefinition.materialAssignments.push_back(assignment);
            assignElement = assignElement->NextSiblingElement("assign");
        }
    }

    const tinyxml2::XMLElement *physicsElement = caseElement->FirstChildElement("physics");
    while (physicsElement != NULL) {
        PhysicsDefinition physics;
        if (physicsElement->Attribute("kind") != NULL) {
            physics.kind = physicsElement->Attribute("kind");
        }
        if (physicsElement->Attribute("order") != NULL) {
            physics.order = std::atoi(physicsElement->Attribute("order"));
            if (physics.order < 1) {
                physics.order = 1;
            }
        }

        // Parse solver configuration
        const tinyxml2::XMLElement *solverElement = physicsElement->FirstChildElement("solver");
        if (solverElement != NULL) {
            if (solverElement->Attribute("type") != NULL) {
                physics.solver.type = solverElement->Attribute("type");
            }
            if (solverElement->Attribute("max_iter") != NULL) {
                physics.solver.maxIterations = std::atoi(solverElement->Attribute("max_iter"));
            }
            if (solverElement->Attribute("tolerance") != NULL) {
                physics.solver.relativeTolerance = std::atof(solverElement->Attribute("tolerance"));
            }
            if (solverElement->Attribute("print_level") != NULL) {
                physics.solver.printLevel = std::atoi(solverElement->Attribute("print_level"));
            }
        }

        const tinyxml2::XMLElement *boundaryElement = physicsElement->FirstChildElement("boundary");
        while (boundaryElement != NULL) {
            BoundaryCondition boundary;
            if (boundaryElement->Attribute("kind") != NULL) {
                boundary.kind = boundaryElement->Attribute("kind");
            }
            if (boundaryElement->Attribute("value") != NULL) {
                boundary.valueText = boundaryElement->Attribute("value");
            }
            if (boundaryElement->Attribute("aux") != NULL) {
                boundary.auxText = boundaryElement->Attribute("aux");
            }
            if (!parseIdAttribute(boundaryElement->Attribute("ids"), boundary.ids, errorMessage)) {
                errorMessage = "Failed to parse boundary ids: " + errorMessage;
                return false;
            }
            physics.boundaries.push_back(boundary);
            boundaryElement = boundaryElement->NextSiblingElement("boundary");
        }

        const tinyxml2::XMLElement *sourceElement = physicsElement->FirstChildElement("source");
        while (sourceElement != NULL) {
            SourceDefinition source;
            if (sourceElement->Attribute("kind") != NULL) {
                source.kind = sourceElement->Attribute("kind");
            }
            if (sourceElement->Attribute("value") != NULL) {
                source.valueText = sourceElement->Attribute("value");
            }
            if (!parseIdAttribute(sourceElement->Attribute("domains"), source.domainIds, errorMessage)) {
                errorMessage = "Failed to parse source domains: " + errorMessage;
                return false;
            }
            physics.sources.push_back(source);
            sourceElement = sourceElement->NextSiblingElement("source");
        }

        caseDefinition.physicsDefinitions.push_back(physics);
        physicsElement = physicsElement->NextSiblingElement("physics");
    }

    const tinyxml2::XMLElement *coupledElement = caseElement->FirstChildElement("coupledPhysics");
    while (coupledElement != NULL) {
        CoupledPhysicsDefinition coupling;
        if (coupledElement->Attribute("name") != NULL) {
            coupling.name = coupledElement->Attribute("name");
        }
        if (coupledElement->Attribute("kind") != NULL) {
            coupling.kind = coupledElement->Attribute("kind");
        }
        if (coupledElement->Attribute("physics") != NULL) {
            splitCsv(coupledElement->Attribute("physics"), coupling.physicsKinds);
        }
        if (!parseIdAttribute(coupledElement->Attribute("domains"), coupling.domainIds, errorMessage)) {
            errorMessage = "Failed to parse coupled physics domains: " + errorMessage;
            return false;
        }
        caseDefinition.coupledPhysicsDefinitions.push_back(coupling);
        coupledElement = coupledElement->NextSiblingElement("coupledPhysics");
    }

    // Parse coupling configuration (new format)
    const tinyxml2::XMLElement *couplingConfigElement = caseElement->FirstChildElement("coupling");
    if (couplingConfigElement != NULL) {
        if (couplingConfigElement->Attribute("method") != NULL) {
            caseDefinition.couplingConfig.method = couplingConfigElement->Attribute("method");
        }
        if (couplingConfigElement->Attribute("max_iter") != NULL) {
            caseDefinition.couplingConfig.maxIterations = std::atoi(couplingConfigElement->Attribute("max_iter"));
        }
        if (couplingConfigElement->Attribute("tolerance") != NULL) {
            caseDefinition.couplingConfig.tolerance = std::atof(couplingConfigElement->Attribute("tolerance"));
        }
    }

    return true;
}

} // namespace mpfem
