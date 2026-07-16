#ifndef IPARAM_PROXY_H
#define IPARAM_PROXY_H

#include "scalar.h"

/**
 * @file IParamAdapter.h
 * @brief Minimal proxy interface to query parameter-like values by (block, id).
 *
 * This header defines @ref IParameterProxy, a light abstraction used by
 * higher-level components to read parameter values without directly depending
 * on the full Parameters subsystem.
 *
 * Typical implementations wrap providers/adapters such as:
 * - @ref ParameterProvider (via @ref ParameterProxy),
 * - or MARTY-specific accessors (via @ref MartyParameterProxy).
 *
 * @tparam T Block identifier type (commonly std::string).
 * @tparam V Entry identifier type (commonly LhaID).
 */

/**
 * @class IParameterProxy
 * @ingroup DataProvidersModule
 * @brief Abstract interface for reading parameters and basic metadata.
 *
 * Defines a minimal contract to:
 * - retrieve a scalar value for (block, id),
 * - check existence of (block, id),
 * - retrieve a block scale factor.
 *
 * This interface is intentionally small so that clients can be agnostic
 * of the underlying storage (Parameters, external libs, etc.).
 *
 * @tparam T Block identifier type.
 * @tparam V Parameter identifier type.
 */
template<typename T, typename V>
class IParameterProxy {
public:
    /**
     * @brief Retrieves the scalar value for (block, id).
     *
     * @param block Block identifier.
     * @param id    Parameter identifier within the block.
     * @return The scalar value.
     */
    virtual scalar_t operator()(const T& x, const V& y) const = 0;

    /**
     * @brief Checks if a parameter exists for (block, id).
     *
     * @param block Block identifier.
     * @param id    Parameter identifier.
     * @return true if the entry exists, false otherwise.
     */
    virtual bool exist(const T& block, const V& id) const = 0;

    /**
     * @brief Returns the scale associated to a block.
     *
     * Some blocks may define a scaling factor used to interpret values.
     *
     * @param block Block identifier.
     * @return Scale factor for that block.
     */
    virtual double get_scale(const T& block) const = 0;
};

#endif // IPARAM_PROXY_H