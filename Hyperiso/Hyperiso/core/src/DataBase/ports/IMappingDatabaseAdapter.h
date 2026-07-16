#ifndef IMAPPINGDATABASEADAPTER_H
#define IMAPPINGDATABASEADAPTER_H

#include <unordered_map>
#include <string>
#include <memory>
#include <optional>

#include "IParamMappingSource.h"

/**
 * @file IMappingDatabaseAdapter.h
 * @brief Interface for read-only access to a mapping database of parameters.
 *
 * This abstraction sits on top of some persistence mechanism
 * (file-backed, in-memory, SQL, etc.) and exposes an API in terms of
 * InterpretedParam objects.
 */

/**
 * @class IMappingDatabaseAdapter
 * @brief Abstract view of a parameter mapping database.
 *
 * A mapping database adapter is responsible for:
 *   - providing the full parameter map (name -> InterpretedParam),
 *   - offering single-parameter lookup by name,
 *   - reporting a human-readable instance name for diagnostics.
 *
 * Concrete implementations:
 *   - DefaultMappingDatabaseAdapter: thin wrapper around MappingDatabase.
 */
class IMappingDatabaseAdapter {
public:
    virtual ~IMappingDatabaseAdapter() = default;

    /**
     * @brief Returns all parameter mappings stored in the database.
     *
     * @return Map from parameter name to InterpretedParam.
     */
    virtual std::unordered_map<std::string, InterpretedParam>
    getParams() const = 0;

    /**
     * @brief Returns a single parameter mapping by name.
     *
     * @param name Logical parameter name (e.g. "m_B_s", "alpha_s").
     * @return Optional InterpretedParam; std::nullopt if the parameter
     *         is not present in the database.
     */
    virtual std::optional<InterpretedParam>
    getParam(const std::string& name) const = 0;

    /**
     * @brief Returns a human-readable name for this mapping instance.
     *
     * This is useful when working with multiple mapping databases at once
     * (e.g. "default", "benchmark1", "user_override", ...).
     *
     * @return Instance name.
     */
    virtual std::string instanceName() const = 0;
};

#endif // IMAPPINGDATABASEADAPTER_H
