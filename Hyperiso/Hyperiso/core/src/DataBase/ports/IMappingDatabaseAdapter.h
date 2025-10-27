#ifndef IMAPPINGDATABASEADAPTER_H
#define IMAPPINGDATABASEADAPTER_H

#include <unordered_map>
#include <string>
#include <memory>
#include <optional>
#include "IParamMappingSource.h"

class IMappingDatabaseAdapter {
public:
    virtual ~IMappingDatabaseAdapter() = default;

    virtual std::unordered_map<std::string, InterpretedParam>
    getParams() const = 0;

    virtual std::optional<InterpretedParam>
    getParam(const std::string& name) const = 0;

    virtual std::string instanceName() const = 0;
};

#endif // IMAPPINGDATABASEADAPTER_H
