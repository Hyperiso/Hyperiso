#include "MappingDatabaseProxy.h"
#include <stdexcept>

std::shared_ptr<MappingDatabaseProxy> MappingDatabaseProxy::fromFactory(
    const IMappingAdapterFactory& factory,
    const std::string& instanceName,
    const std::string& jsonFilePath,
    std::shared_ptr<IParamMappingSource> loader
) {
    auto adapter = factory.create(instanceName, jsonFilePath, std::move(loader));
    return std::make_shared<MappingDatabaseProxy>(std::move(adapter));
}

MappingDatabaseProxy::MappingDatabaseProxy(std::shared_ptr<IMappingDatabaseAdapter> adapter)
    : adapter_(std::move(adapter)) {
    if (!adapter_) throw std::runtime_error("MappingDatabaseProxy: adapter nul");
}

const std::unordered_map<std::string, InterpretedParam>&
MappingDatabaseProxy::getParams() const {
    return adapter_->getParams();
}

std::optional<InterpretedParam>
MappingDatabaseProxy::getParam(const std::string& name) const {
    return adapter_->getParam(name);
}

std::string MappingDatabaseProxy::instanceName() const {
    return adapter_->instanceName();
}
