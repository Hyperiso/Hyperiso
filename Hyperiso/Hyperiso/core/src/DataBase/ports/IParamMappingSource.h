#ifndef IPARAMMAPPINGSOURCE_H
#define IPARAMMAPPINGSOURCE_H

#include <memory>
#include <string>
#include <unordered_map>
#include "Include.h"
#include "IParser.h"
#include "InterpretedParam.h"

class IParamMappingSource {
public:
    virtual ~IParamMappingSource() = default;

    virtual std::unordered_map<std::string, InterpretedParam>
    loadFromFile(const std::string& filePath) const = 0;
};

#endif // IPARAMMAPPINGSOURCE_H
