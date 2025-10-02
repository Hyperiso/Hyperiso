#ifndef DEFAULTMAPPINGADAPTERFACTORY_H
#define DEFAULTMAPPINGADAPTERFACTORY_H

#include "IMappingAdapterFactory.h"
#include "DefaultMappingDatabaseAdapter.h"
#include "JsonParamMappingAdapter.h"
#include "JsonParser.h"
#include "MappingDatabase.h"

class DefaultMappingAdapterFactory : public IMappingAdapterFactory {
public:
    std::shared_ptr<IMappingDatabaseAdapter> create(
        const std::string& instanceName,
        const std::string& jsonFilePath,
        std::shared_ptr<IParamMappingSource> loader
    ) const override
    {
        if (!loader) {
            loader = std::make_shared<JsonParamMappingAdapter>(std::make_shared<JSONParser>());
        }

        auto db = MappingDatabase::getInstance(instanceName, jsonFilePath, loader);
        if (!db) throw std::runtime_error("DefaultMappingAdapterFactory: DB introuvable");

        return std::make_shared<DefaultMappingDatabaseAdapter>(instanceName, db);
    }
};

#endif
