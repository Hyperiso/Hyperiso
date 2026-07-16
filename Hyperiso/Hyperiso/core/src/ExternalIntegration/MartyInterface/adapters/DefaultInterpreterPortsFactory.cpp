#include "DefaultInterpreterPortsFactory.h"

#include "MappingDatabaseProxy.h"
#include "ParameterResolverProxy.h"
#include "DefaultMappingAdapterFactory.h"
#include "JsonParamMappingAdapter.h"
#include "JsonParser.h"

#include <filesystem>
#include <sstream>
#include <stdexcept>
#include <system_error>

namespace {

namespace fs = std::filesystem;

bool same_mapping_file(const std::string& lhs, const std::string& rhs)
{
    if (lhs == rhs) {
        return true;
    }

    std::error_code ec;
    const fs::path lhs_path(lhs);
    const fs::path rhs_path(rhs);

    if (fs::exists(lhs_path, ec) && fs::exists(rhs_path, ec)) {
        ec.clear();
        if (fs::equivalent(lhs_path, rhs_path, ec) && !ec) {
            return true;
        }
    }

    return fs::absolute(lhs_path).lexically_normal() == fs::absolute(rhs_path).lexically_normal();
}

void validate_no_mapping_collision(const IMappingDatabasePort& smDB,
                                   const IMappingDatabasePort& modelDB,
                                   const std::string& modelName,
                                   const std::string& modelJsonPath,
                                   const std::string& smJsonPath)
{
    if (same_mapping_file(modelJsonPath, smJsonPath)) {
        return;
    }

    const auto smParams = smDB.getParams();
    const auto modelParams = modelDB.getParams();

    for (const auto& [name, _] : modelParams) {
        if (smParams.contains(name)) {
            std::ostringstream oss;
            oss << "MARTY mapping collision for key '" << name << "' while loading model '"
                << modelName << "'. The key exists in both the read-only SM mapping '"
                << smJsonPath << "' and the user/model mapping '" << modelJsonPath
                << "'. Rename the BSM key or remove the duplicate entry.";
            throw std::runtime_error(oss.str());
        }
    }
}

std::string resolve_sm_mapping_path(const std::string& modelName,
                                    const std::string& modelJsonPath,
                                    const std::string& smJsonPath)
{
    if (modelName == "SM" || !same_mapping_file(modelJsonPath, smJsonPath)) {
        return smJsonPath;
    }

    // Backward-compatibility guard: older call sites may still pass
    // FileNameManager::getjsondbmodel() for both the model mapping and the SM
    // mapping. When a BSM mapping is configured, that makes both databases load
    // the BSM JSON and SM parameters such as G_F become unavailable. In that
    // case, recover the read-only SM mapping from the same mapping directory.
    std::error_code ec;
    const fs::path candidate = fs::path(modelJsonPath).parent_path() / "sm.json";
    if (fs::exists(candidate, ec) && !ec && !same_mapping_file(candidate.string(), modelJsonPath)) {
        return candidate.string();
    }

    return smJsonPath;
}

} // namespace

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

    const std::string resolvedSmJsonPath = resolve_sm_mapping_path(modelName, modelJsonPath, smJsonPath);

    auto modelDB = MappingDatabaseProxy::fromFactory(*adapter, modelName, modelJsonPath, src);
    auto smDB    = MappingDatabaseProxy::fromFactory(*adapter, "SM",    resolvedSmJsonPath,    src);

    validate_no_mapping_collision(*smDB, *modelDB, modelName, modelJsonPath, resolvedSmJsonPath);

    return std::make_unique<ParameterResolverProxy>(smDB, modelDB);
}
