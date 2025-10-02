#include "DefaultInterpreterPortsFactory.h"

#include "MappingDatabaseProxy.h"
#include "ParameterResolverProxy.h"
#include "DefaultMappingAdapterFactory.h"
#include "JsonParamMappingAdapter.h"
#include "JsonParser.h"

DefaultInterpreterPortsFactory::DefaultInterpreterPortsFactory(
    std::shared_ptr<IMappingAdapterFactory> adapterFactory,
    std::shared_ptr<IParamMappingSource>    loader
)
: adapterFactory_(std::move(adapterFactory))
, loader_(std::move(loader))
{}

std::unique_ptr<IParameterResolver>
DefaultInterpreterPortsFactory::makeResolver(const std::string& modelName,
                                             const std::string& modelJsonPath,
                                             const std::string& smJsonPath) const
{
    std::shared_ptr<IMappingAdapterFactory> adapter =
        adapterFactory_ ? adapterFactory_
                        : std::make_shared<DefaultMappingAdapterFactory>();

    std::shared_ptr<IParamMappingSource> src =
        loader_ ? loader_
                : std::make_shared<JsonParamMappingAdapter>(std::make_shared<JSONParser>());

    auto modelDB = MappingDatabaseProxy::fromFactory(*adapter, modelName, modelJsonPath, src);
    auto smDB    = MappingDatabaseProxy::fromFactory(*adapter, "SM",    smJsonPath,    src);

    return std::make_unique<ParameterResolverProxy>(smDB, modelDB);
}
