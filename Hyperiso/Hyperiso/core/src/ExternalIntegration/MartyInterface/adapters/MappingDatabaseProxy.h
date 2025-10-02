#ifndef MAPPING_DATABASE_PROXY_H
#define MAPPING_DATABASE_PROXY_H

#include <memory>
#include <string>
#include <unordered_map>
#include <optional>

#include "IMappingDatabasePort.h"
#include "IMappingDatabaseAdapter.h"
#include "IMappingAdapterFactory.h"

class MappingDatabaseProxy final : public IMappingDatabasePort {
public:
    static std::shared_ptr<MappingDatabaseProxy> fromFactory(
        const IMappingAdapterFactory& factory,
        const std::string& instanceName,
        const std::string& jsonFilePath,
        std::shared_ptr<IParamMappingSource> loader = nullptr
    );

    explicit MappingDatabaseProxy(std::shared_ptr<IMappingDatabaseAdapter> adapter);

    const std::unordered_map<std::string, InterpretedParam>& getParams() const override;
    std::optional<InterpretedParam> getParam(const std::string& name) const override;
    std::string instanceName() const override;

private:
    std::shared_ptr<IMappingDatabaseAdapter> adapter_;
};

#endif // MAPPING_DATABASE_PROXY_H
