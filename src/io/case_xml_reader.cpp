#include "case_xml_reader.hpp"
#include "boundary_spec.hpp"
#include "value_parser.hpp"
#include "logger.hpp"

#include "general/tinyxml2.h"

#include <sstream>
#include <string>
#include <cstdlib>

namespace mpfem
{

    namespace
    {
        const char TOKEN_DELIMITER = ',';
        const char RANGE_DELIMITER = '-';
        const int MIN_VALID_ID = 1;

        std::string trim(const std::string &text)
        {
            std::size_t first = 0;
            while (first < text.size() && std::isspace(static_cast<unsigned char>(text[first])) != 0)
            {
                ++first;
            }

            std::size_t last = text.size();
            while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1])) != 0)
            {
                --last;
            }

            return text.substr(first, last - first);
        }

        bool parsePositiveInt(const std::string &text, int &value)
        {
            if (text.empty())
            {
                return false;
            }

            for (std::size_t i = 0; i < text.size(); ++i)
            {
                if (std::isdigit(static_cast<unsigned char>(text[i])) == 0)
                {
                    return false;
                }
            }

            try
            {
                value = std::stoi(text);
            }
            catch (...)
            {
                return false;
            }

            return value >= MIN_VALID_ID;
        }

        void parseIds(const std::string &text, std::set<int> &ids)
        {
            ids.clear();

            std::stringstream tokenStream(text);
            std::string token;

            while (std::getline(tokenStream, token, TOKEN_DELIMITER))
            {
                const std::string trimmedToken = trim(token);
                Check(!trimmedToken.empty(), "Empty token found in range expression");

                const std::size_t rangePos = trimmedToken.find(RANGE_DELIMITER);
                if (rangePos == std::string::npos)
                {
                    int value = 0;
                    Check(parsePositiveInt(trimmedToken, value),
                          ("Invalid identifier token: " + trimmedToken).c_str());
                    ids.insert(value);
                    continue;
                }

                const std::string beginToken = trim(trimmedToken.substr(0, rangePos));
                const std::string endToken = trim(trimmedToken.substr(rangePos + 1));
                int beginValue = 0;
                int endValue = 0;
                Check(parsePositiveInt(beginToken, beginValue) && parsePositiveInt(endToken, endValue),
                      ("Invalid range token: " + trimmedToken).c_str());
                Check(beginValue <= endValue,
                      ("Descending range is not allowed: " + trimmedToken).c_str());

                for (int value = beginValue; value <= endValue; ++value)
                {
                    ids.insert(value);
                }
            }

            Check(!ids.empty(), "No ids found in expression");
        }

        void parseIdAttribute(const char *attributeText, std::set<int> &target)
        {
            if (attributeText == nullptr)
            {
                target.clear();
                return;
            }
            parseIds(attributeText, target);
        }

        void splitCsv(const std::string &csv, std::vector<std::string> &tokens)
        {
            tokens.clear();
            std::stringstream stream(csv);
            std::string item;
            while (std::getline(stream, item, ','))
            {
                std::string trimmed = item;
                while (!trimmed.empty() && (trimmed.front() == ' ' || trimmed.front() == '\t'))
                {
                    trimmed.erase(trimmed.begin());
                }
                while (!trimmed.empty() && (trimmed.back() == ' ' || trimmed.back() == '\t'))
                {
                    trimmed.pop_back();
                }
                if (!trimmed.empty())
                {
                    tokens.push_back(trimmed);
                }
            }
        }

    } // namespace

    void CaseXmlReader::readFromFile(const std::string &filePath, CaseDefinition &caseDefinition)
    {
        caseDefinition = CaseDefinition();

        tinyxml2::XMLDocument document;
        const tinyxml2::XMLError loadError = document.LoadFile(filePath.c_str());
        Check(loadError == tinyxml2::XML_SUCCESS, "Failed to parse XML file: " + filePath);

        const tinyxml2::XMLElement *caseElement = document.FirstChildElement("case");
        Check(caseElement != nullptr, "Missing <case> root node");

        if (caseElement->Attribute("name") != nullptr)
        {
            caseDefinition.caseName = caseElement->Attribute("name");
        }

        const tinyxml2::XMLElement *studyElement = caseElement->FirstChildElement("study");
        if (studyElement != nullptr && studyElement->Attribute("type") != nullptr)
        {
            caseDefinition.studyType = studyElement->Attribute("type");
        }

        const tinyxml2::XMLElement *pathsElement = caseElement->FirstChildElement("paths");
        if (pathsElement != nullptr)
        {
            if (pathsElement->Attribute("mesh") != nullptr)
            {
                caseDefinition.meshPath = pathsElement->Attribute("mesh");
            }
            if (pathsElement->Attribute("materials") != nullptr)
            {
                caseDefinition.materialsPath = pathsElement->Attribute("materials");
            }
            if (pathsElement->Attribute("comsol_result") != nullptr)
            {
                caseDefinition.comsolResultPath = pathsElement->Attribute("comsol_result");
            }
        }

        const tinyxml2::XMLElement *variablesElement = caseElement->FirstChildElement("variables");
        if (variablesElement != nullptr)
        {
            const tinyxml2::XMLElement *variableElement = variablesElement->FirstChildElement("var");
            while (variableElement != nullptr)
            {
                VariableEntry entry;
                if (variableElement->Attribute("name") != nullptr)
                {
                    entry.name = variableElement->Attribute("name");
                }
                if (variableElement->Attribute("value") != nullptr)
                {
                    entry.valueText = variableElement->Attribute("value");
                }
                if (variableElement->Attribute("si") != nullptr)
                {
                    ValueParser::parseFirstNumber(variableElement->Attribute("si"), entry.siValue);
                }
                caseDefinition.variables.push_back(entry);
                variableElement = variableElement->NextSiblingElement("var");
            }
        }

        const tinyxml2::XMLElement *materialsElement = caseElement->FirstChildElement("materials");
        if (materialsElement != nullptr)
        {
            const tinyxml2::XMLElement *assignElement = materialsElement->FirstChildElement("assign");
            while (assignElement != nullptr)
            {
                MaterialAssignment assignment;
                if (assignElement->Attribute("material") != nullptr)
                {
                    assignment.materialTag = assignElement->Attribute("material");
                }
                parseIdAttribute(assignElement->Attribute("domains"), assignment.domainIds);
                caseDefinition.materialAssignments.push_back(assignment);
                assignElement = assignElement->NextSiblingElement("assign");
            }
        }

        const tinyxml2::XMLElement *physicsElement = caseElement->FirstChildElement("physics");
        while (physicsElement != nullptr)
        {
            PhysicsDefinition physics;
            if (physicsElement->Attribute("kind") != nullptr)
            {
                physics.kind = physicsElement->Attribute("kind");
            }
            if (physicsElement->Attribute("order") != nullptr)
            {
                physics.order = std::atoi(physicsElement->Attribute("order"));
                if (physics.order < 1)
                {
                    physics.order = 1;
                }
            }

            // Parse solver configuration
            const tinyxml2::XMLElement *solverElement = physicsElement->FirstChildElement("solver");
            if (solverElement != nullptr)
            {
                if (solverElement->Attribute("type") != nullptr)
                {
                    physics.solver.type = solverElement->Attribute("type");
                }
                if (solverElement->Attribute("max_iter") != nullptr)
                {
                    physics.solver.maxIterations = std::atoi(solverElement->Attribute("max_iter"));
                }
                if (solverElement->Attribute("tolerance") != nullptr)
                {
                    physics.solver.relativeTolerance = std::atof(solverElement->Attribute("tolerance"));
                }
                if (solverElement->Attribute("print_level") != nullptr)
                {
                    physics.solver.printLevel = std::atoi(solverElement->Attribute("print_level"));
                }
            }

            const tinyxml2::XMLElement *boundaryElement = physicsElement->FirstChildElement("boundary");
            while (boundaryElement != nullptr)
            {
                BoundaryCondition boundary;
                if (boundaryElement->Attribute("kind") != nullptr)
                {
                    boundary.kind = boundaryElement->Attribute("kind");
                }
                parseIdAttribute(boundaryElement->Attribute("ids"), boundary.ids);

                // Parse parameters from child <param> elements
                const tinyxml2::XMLElement *paramElement = boundaryElement->FirstChildElement("param");
                while (paramElement != nullptr)
                {
                    const char *nameAttr = paramElement->Attribute("name");
                    const char *valueAttr = paramElement->Attribute("value");
                    if (nameAttr != nullptr && valueAttr != nullptr)
                    {
                        boundary.params[nameAttr] = valueAttr;
                    }
                    paramElement = paramElement->NextSiblingElement("param");
                }

                physics.boundaries.push_back(boundary);
                boundaryElement = boundaryElement->NextSiblingElement("boundary");
            }

            const tinyxml2::XMLElement *sourceElement = physicsElement->FirstChildElement("source");
            while (sourceElement != nullptr)
            {
                SourceDefinition source;
                if (sourceElement->Attribute("kind") != nullptr)
                {
                    source.kind = sourceElement->Attribute("kind");
                }
                if (sourceElement->Attribute("value") != nullptr)
                {
                    source.valueText = sourceElement->Attribute("value");
                }
                parseIdAttribute(sourceElement->Attribute("domains"), source.domainIds);
                physics.sources.push_back(source);
                sourceElement = sourceElement->NextSiblingElement("source");
            }

            caseDefinition.physicsDefinitions.push_back(physics);
            physicsElement = physicsElement->NextSiblingElement("physics");
        }

        const tinyxml2::XMLElement *coupledElement = caseElement->FirstChildElement("coupledPhysics");
        while (coupledElement != nullptr)
        {
            CoupledPhysicsDefinition coupling;
            if (coupledElement->Attribute("name") != nullptr)
            {
                coupling.name = coupledElement->Attribute("name");
            }
            if (coupledElement->Attribute("kind") != nullptr)
            {
                coupling.kind = coupledElement->Attribute("kind");
            }
            if (coupledElement->Attribute("physics") != nullptr)
            {
                splitCsv(coupledElement->Attribute("physics"), coupling.physicsKinds);
            }
            parseIdAttribute(coupledElement->Attribute("domains"), coupling.domainIds);
            caseDefinition.coupledPhysicsDefinitions.push_back(coupling);
            coupledElement = coupledElement->NextSiblingElement("coupledPhysics");
        }

        // Parse coupling configuration
        const tinyxml2::XMLElement *couplingConfigElement = caseElement->FirstChildElement("coupling");
        if (couplingConfigElement != nullptr)
        {
            if (couplingConfigElement->Attribute("method") != nullptr)
            {
                caseDefinition.couplingConfig.method = couplingConfigElement->Attribute("method");
            }
            if (couplingConfigElement->Attribute("max_iter") != nullptr)
            {
                caseDefinition.couplingConfig.maxIterations = std::atoi(couplingConfigElement->Attribute("max_iter"));
            }
            if (couplingConfigElement->Attribute("tolerance") != nullptr)
            {
                caseDefinition.couplingConfig.tolerance = std::atof(couplingConfigElement->Attribute("tolerance"));
            }
        }
    }

} // namespace mpfem