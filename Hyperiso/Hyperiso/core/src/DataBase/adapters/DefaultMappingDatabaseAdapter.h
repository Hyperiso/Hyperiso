#ifndef DEFAULTMAPPINGDATABASEADAPTER_H
#define DEFAULTMAPPINGDATABASEADAPTER_H

#include <utility>

#include "IMappingDatabaseAdapter.h"
#include "MappingDatabase.h"

/**
 * @file DefaultMappingDatabaseAdapter.h
 * @brief Thin adapter exposing MappingDatabase through IMappingDatabaseAdapter.
 *
 * This adapter simply forwards calls to an underlying MappingDatabase while
 * providing an instance name for identification.
 */

/**
 * @class DefaultMappingDatabaseAdapter
 * @brief Default implementation of IMappingDatabaseAdapter.
 *
 * It wraps a MappingDatabase and:
 *   - delegates getParams() directly,
 *   - implements getParam(name) as a lookup in the full map,
 *   - exposes an instanceName() string for higher-level identification.
 */
class DefaultMappingDatabaseAdapter : public IMappingDatabaseAdapter {
public:
    /**
     * @brief Constructs an adapter around a MappingDatabase instance.
     *
     * @param instanceName Human-readable name for this mapping instance.
     * @param db           Shared pointer to the underlying MappingDatabase.
     */
    explicit DefaultMappingDatabaseAdapter(std::string instanceName,
                                           std::shared_ptr<MappingDatabase> db)
        : name_(std::move(instanceName)), db_(std::move(db)) {}
    
    /**
     * @brief Returns the full parameter map from the backing database.
     */
    std::unordered_map<std::string, InterpretedParam>
    getParams() const override {
        return db_->getParams();
    }

    /**
     * @brief Returns a single parameter mapping by name, if present.
     *
     * @param name Parameter name.
     * @return Optional InterpretedParam, std::nullopt if not found.
     */
    std::optional<InterpretedParam>
    getParam(const std::string& name) const override {
        auto m = getParams();
        if (auto it = m.find(name); it != m.end()) return it->second;
        return std::nullopt;
    }

    /**
     * @brief Returns the logical instance name of this adapter.
     */
    std::string instanceName() const override { return name_; }

private:
    std::string name_;                      ///< Identifier for this instance.
    std::shared_ptr<MappingDatabase> db_;   ///< Underlying persistence object.
};

#endif
