#ifndef MAPPING_DATABASE_PROXY_H
#define MAPPING_DATABASE_PROXY_H

#include <memory>
#include <string>
#include <unordered_map>
#include <optional>

#include "IMappingDatabasePort.h"
#include "IMappingDatabaseAdapter.h"
#include "IMappingAdapterFactory.h"

/**
 * @file MappingDatabaseProxy.h
 * @brief Declares a proxy implementation of IMappingDatabasePort.
 *
 * The ::MappingDatabaseProxy wraps an ::IMappingDatabaseAdapter and exposes
 * it through the ::IMappingDatabasePort interface, providing a stable view
 * to clients regardless of the underlying data source.
 */

/**
 * @class MappingDatabaseProxy
 * @ingroup MappingModule
 * @brief Concrete IMappingDatabasePort that delegates to an adapter.
 *
 * This class:
 *  - is created via ::MappingDatabaseProxy::fromFactory, which
 *    uses an ::IMappingAdapterFactory to build a concrete adapter,
 *  - forwards all queries to the underlying ::IMappingDatabaseAdapter.
 */
class MappingDatabaseProxy final : public IMappingDatabasePort {
public:
    /**
     * @brief Factory helper to create a proxy with a given adapter factory.
     *
     * @param factory      Adapter factory that knows how to read the mapping file.
     * @param instanceName Logical name of this database instance (e.g. model name).
     * @param jsonFilePath Path to the mapping JSON file.
     * @param loader       Optional mapping source loader (may be null).
     * @return Shared pointer to the newly constructed proxy.
     */
    static std::shared_ptr<MappingDatabaseProxy> fromFactory(
        const IMappingAdapterFactory& factory,
        const std::string& instanceName,
        const std::string& jsonFilePath,
        std::shared_ptr<IParamMappingSource> loader = nullptr
    );

    /**
     * @brief Constructs a new proxy from an existing adapter.
     * @param adapter Concrete mapping database adapter (must be non-null).
     *
     * @throws std::runtime_error If @p adapter is null.
     */
    explicit MappingDatabaseProxy(std::shared_ptr<IMappingDatabaseAdapter> adapter);

    /// @copydoc IMappingDatabasePort::getParams()
    std::unordered_map<std::string, InterpretedParam> getParams() const override;

    /// @copydoc IMappingDatabasePort::getParam()
    std::optional<InterpretedParam> getParam(const std::string& name) const override;

    /// @copydoc IMappingDatabasePort::instanceName()
    std::string instanceName() const override;

private:
    /// Underlying adapter that actually stores and parses the mapping.
    std::shared_ptr<IMappingDatabaseAdapter> adapter_;
};

#endif // MAPPING_DATABASE_PROXY_H
