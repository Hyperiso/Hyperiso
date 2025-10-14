#ifndef IMAPPING_DATABASE_PORT_H
#define IMAPPING_DATABASE_PORT_H

#include <memory>
#include <string>
#include <unordered_map>
#include <optional>
#include "IParamMappingSource.h"
#include "InterpretedParam.h"

class IMappingDatabasePort {
public:
    virtual ~IMappingDatabasePort() = default;

    virtual std::unordered_map<std::string, InterpretedParam>
    getParams() const = 0;

    virtual std::optional<InterpretedParam>
    getParam(const std::string& name) const = 0;

    virtual std::string instanceName() const = 0;
};

#endif