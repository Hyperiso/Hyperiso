// IParamMappingSource.h
#ifndef IPARAMMAPPINGSOURCE_H
#define IPARAMMAPPINGSOURCE_H

#include <memory>
#include <string>
#include <unordered_map>
#include "Include.h" // pour LhaID
#include "IParser.h" // interface de ton parser

struct InterpretedParam {
    std::string block;
    LhaID pdgCode;
};

class IParamMappingSource {
public:
    virtual ~IParamMappingSource() = default;

    // Charge depuis un fichier JSON (ou autre), retourne la map prête à l’emploi.
    virtual std::unordered_map<std::string, InterpretedParam>
    loadFromFile(const std::string& filePath) const = 0;
};

#endif // IPARAMMAPPINGSOURCE_H
