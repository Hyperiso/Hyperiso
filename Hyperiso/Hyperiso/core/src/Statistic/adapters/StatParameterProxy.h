#ifndef STAT_PARAMETER_PROXY_H
#define STAT_PARAMETER_PROXY_H

#include <memory>
#include <string>
#include <unordered_set>

#include "IStatParameterProxy.h"
#include "ParameterProvider.h"
#include "Include.h"
#include "BinnedObservableId.h"

/**
 * @file StatParameterProxy.h
 * @brief Statistics-layer proxy for read-only access to parameters and observables.
 *
 * This proxy adapts @ref ParameterProvider to the needs of the statistics module.
 * It provides:
 * - access to generic typed parameters via @ref ParamId,
 * - access to typed block entries via `(block, LhaID)`,
 * - access to observable entries via @ref ObservableId,
 * - access to binned observable entries via @ref BinnedObservableId,
 * - and direct retrieval of the underlying @ref Parameter objects.
 *
 * Two internal providers are maintained:
 * - `pp`           : generic provider used with fully typed @ref ParamId access,
 * - `pp_with_type` : provider bound to a specific @ref ParameterType, used for
 *                    `(block, LhaID)` queries.
 *
 * Special handling is provided for Wilson parameters:
 * if a Wilson entry is missing, the proxy returns a default-constructed
 * @ref scalar_t instead of failing. This is convenient because some Wilson
 * blocks may be optional, backend-dependent, or group-dependent.
 *
 * @see IStatParameterProxy
 * @see ParameterProvider
 * @see Parameter
 */

/**
 * @class StatParameterProxy
 * @brief Concrete read-only proxy used by the statistics layer.
 *
 * The proxy enforces a whitelist of parameter types accessible from the
 * statistics/business layer:
 * - SM
 * - BSM
 * - WILSON
 * - FLAVOR
 * - DECAY
 * - OBSERVABLE
 *
 * It can be constructed either:
 * - untyped in practice through the generic provider, while still exposing
 *   typed @ref ParamId access,
 * - or bound to a given @ref ParameterType for block-based queries.
 */
class StatParameterProxy : public IStatParameterProxy {
public:
    /**
     * @brief Constructs a proxy optionally bound to a parameter type.
     *
     * If @p type is not allowed by the statistics layer, a ValueError is logged.
     *
     * @param type Parameter namespace handled by block-based access.
     */
    StatParameterProxy(ParameterType type = ParameterType::SM);

    /**
     * @brief Returns the parameter object identified by a typed @ref ParamId.
     *
     * @param pid Fully typed parameter identifier.
     * @return Shared pointer to the corresponding @ref Parameter.
     */
    std::shared_ptr<Parameter> get_param(const ParamId&) const override;

    /**
     * @brief Returns the parameter object identified by block and LHA id.
     *
     * This uses the internal provider bound to the constructor-selected
     * @ref ParameterType.
     *
     * @param block Block name.
     * @param id    LHA-like identifier.
     * @return Shared pointer to the corresponding @ref Parameter.
     */
    std::shared_ptr<Parameter> get_param(const std::string& block, const LhaID& id) const override;

    /**
     * @brief Returns a scalar quantity for a typed parameter.
     *
     * If the parameter has no type, a warning is emitted.
     * If the type is not allowed, a ValueError is logged.
     *
     * For Wilson parameters, missing entries return a default-constructed
     * @ref scalar_t.
     *
     * @param pid    Fully typed parameter identifier.
     * @param d_type Requested quantity type.
     * @return Requested scalar quantity.
     */
    scalar_t operator()(const ParamId&, DataType d_type=DataType::VALUE) const override;

    /**
     * @brief Returns a scalar quantity for an observable entry.
     *
     * The observable is internally mapped to the observable block (`FOBS`)
     * using its FLHA identifier.
     *
     * @param id     Observable identifier.
     * @param d_type Requested quantity type.
     * @return Requested scalar quantity.
     */
    double operator()(const ObservableId&, DataType d_type=DataType::VALUE) const override;

    /**
     * @brief Returns a scalar quantity for a binned observable entry.
     *
     * The binned observable is internally mapped to the observable block (`FOBS`)
     * using its FLHA identifier.
     *
     * @param id     Binned observable identifier.
     * @param d_type Requested quantity type.
     * @return Requested scalar quantity.
     */
    double operator()(const BinnedObservableId&, DataType d_type=DataType::VALUE) const override;

    /**
     * @brief Returns a scalar quantity for a parameter identified by block and LHA id.
     *
     * This overload uses the provider bound to the constructor-selected
     * @ref ParameterType.
     *
     * For Wilson parameters, missing entries return a default-constructed
     * @ref scalar_t.
     *
     * @param block  Block name.
     * @param id     LHA-like identifier.
     * @param d_type Requested quantity type.
     * @return Requested scalar quantity.
     */
    scalar_t operator()(const std::string& block, const LhaID& id, DataType d_type=DataType::VALUE) const override;

    /**
     * @brief Returns the parameter object associated with a binned observable.
     *
     * The binned observable is resolved to the `FOBS` observable block.
     *
     * @param id Binned observable identifier.
     * @return Shared pointer to the corresponding observable @ref Parameter.
     */
    std::shared_ptr<Parameter> get_obs_param(const BinnedObservableId&) const override;
private:
    ParameterProvider pp;                                               /// Generic provider used for fully typed ParamId access.
    ParameterProvider pp_with_type;                                     /// Provider bound to a specific ParameterType, used for (block, LhaID) access.
    static inline const std::unordered_set<ParameterType> ALLOWED       /// Parameter types allowed through the statistics-layer proxy.
    {ParameterType::SM, ParameterType::BSM, ParameterType::WILSON, ParameterType::FLAVOR, ParameterType::DECAY, ParameterType::OBSERVABLE};
};

#endif 