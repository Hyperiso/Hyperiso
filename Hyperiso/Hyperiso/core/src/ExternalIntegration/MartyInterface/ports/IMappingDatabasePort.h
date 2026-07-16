#ifndef IMAPPING_DATABASE_PORT_H
#define IMAPPING_DATABASE_PORT_H

#include <memory>
#include <string>
#include <unordered_map>
#include <optional>

#include "IParamMappingSource.h"
#include "InterpretedParam.h"

/**
 * @file IMappingDatabasePort.h
 * @brief Declares an abstract port for parameter-mapping databases.
 *
 * This interface abstracts access to a parameter mapping database,
 * decoupling the rest of the system from the concrete storage format
 * (JSON, database, etc.).
 */

/**
 * @class IMappingDatabasePort
 * @ingroup MappingModule
 * @brief Interface for accessing a parameter-mapping database.
 *
 * The interface provides:
 *  - bulk retrieval of all mapped parameters,
 *  - lookup by symbolic name,
 *  - access to a human-readable instance name (e.g. `"SM"`, `"THDM"`).
 */
class IMappingDatabasePort {
public:
    virtual ~IMappingDatabasePort() = default;

    /**
     * @brief Returns all parameters stored in the database.
     *
     * @return Map from parameter name to ::InterpretedParam.
     */
    virtual std::unordered_map<std::string, InterpretedParam>
    getParams() const = 0;

    /**
     * @brief Retrieves a parameter mapping by name.
     *
     * @param name Symbolic parameter name as used in the generated code.
     * @return The corresponding ::InterpretedParam if found, `std::nullopt` otherwise.
     */
    virtual std::optional<InterpretedParam>
    getParam(const std::string& name) const = 0;

    /**
     * @brief Returns a logical name for this database instance.
     *
     * Typically this is the model name, e.g. `"SM"` or `"THDM"`.
     *
     * @return Instance name as a string.
     */
    virtual std::string instanceName() const = 0;
};

#endif