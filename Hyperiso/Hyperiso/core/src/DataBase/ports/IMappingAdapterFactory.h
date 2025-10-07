// IMappingAdapterFactory.h
#ifndef IMAPPINGADAPTERFACTORY_H
#define IMAPPINGADAPTERFACTORY_H

#include <memory>
#include <string>
#include "IMappingDatabaseAdapter.h"
#include "IParamMappingSource.h"

class IMappingAdapterFactory {
public:
    virtual ~IMappingAdapterFactory() = default;

    virtual std::shared_ptr<IMappingDatabaseAdapter> create(
        const std::string& instanceName,
        const std::string& jsonFilePath,
        std::shared_ptr<IParamMappingSource> loader /* peut être nul, fallback interne */
    ) const = 0;
};

#endif // IMAPPINGADAPTERFACTORY_H
