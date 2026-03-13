#ifndef ISTAT_PARAMETER_PROXY_H
#define ISTAT_PARAMETER_PROXY_H

#include <memory>
#include <string>

#include "CorrelationProvider.h"
#include "Include.h"

/**
 * @file IStatParameterProxy.h
 * @brief Statistical-layer port for read-only access to parameters and observable entries.
 *
 * This interface defines the read API used by the statistics module to access:
 * - generic model parameters through @ref ParamId,
 * - typed parameters through `(block, LhaID)`,
 * - observable entries through @ref ObservableId,
 * - binned observable entries through @ref BinnedObservableId,
 * - and the underlying @ref Parameter objects when direct metadata access is needed.
 *
 * The proxy is intentionally read-only: it is meant for likelihood construction,
 * covariance handling, pulls / nuisance inspection, and statistical post-processing,
 * without exposing mutation logic.
 *
 * The returned quantity depends on @ref DataType:
 * - @ref DataType::VALUE        : central value
 * - @ref DataType::STD_STAT     : statistical uncertainty
 * - @ref DataType::STD_SYST     : systematic uncertainty
 * - @ref DataType::STD_COMBINED : quadratic combination of uncertainties
 *
 * Concrete implementations are typically backed by @ref ParameterProvider and
 * the @ref Parameters singleton infrastructure.
 *
 * @see StatParameterProxy
 * @see ParameterProvider
 * @see Parameter
 */

/**
 * @class IStatParameterProxy
 * @brief Abstract read-only interface to parameters and observables for the statistics layer.
 *
 * This interface provides two families of access:
 * - value access through overloaded `operator()`,
 * - object access through `get_param()` / `get_obs_param()`.
 *
 * The object-returning methods are useful when the caller needs more than the
 * scalar value itself (uncertainties, binning, metadata, ownership, etc.).
 */
class IStatParameterProxy {
public:
    /**
     * @brief Virtual destructor.
     */
    virtual ~IStatParameterProxy() = default;

    /**
     * @brief Returns the underlying parameter object identified by a typed @ref ParamId.
     *
     * @param pid Fully typed parameter identifier.
     * @return Shared pointer to the corresponding @ref Parameter.
     */
    virtual std::shared_ptr<Parameter> get_param(const ParamId&) const = 0;

    /**
     * @brief Returns the underlying parameter object identified by block and LHA id.
     *
     * This overload is intended for proxies already bound to a specific
     * @ref ParameterType.
     *
     * @param block Block name.
     * @param id    LHA-like identifier inside the block.
     * @return Shared pointer to the corresponding @ref Parameter.
     */
    virtual std::shared_ptr<Parameter> get_param(const std::string& block, const LhaID& id) const = 0;

    /**
     * @brief Returns a scalar quantity for a typed parameter.
     *
     * Depending on @p d_type, returns the central value or one of its uncertainties.
     *
     * @param pid    Fully typed parameter identifier.
     * @param d_type Requested quantity type.
     * @return Requested scalar quantity.
     */
    virtual scalar_t operator()(const ParamId&, DataType d_type=DataType::VALUE) const = 0;

    /**
     * @brief Returns a scalar quantity for an observable identified by @ref ObservableId.
     *
     * This is typically used for experimental observable entries stored in
     * the observable parameter blocks.
     *
     * @param id     Observable identifier.
     * @param d_type Requested quantity type.
     * @return Requested scalar quantity.
     */
    virtual double operator()(const ObservableId&, DataType d_type=DataType::VALUE) const = 0;

    /**
     * @brief Returns a scalar quantity for a binned observable entry.
     *
     * This overload is used when the observable is identified not only by its
     * base observable id but also by a bin definition encoded in
     * @ref BinnedObservableId.
     *
     * @param id     Binned observable identifier.
     * @param d_type Requested quantity type.
     * @return Requested scalar quantity.
     */
    virtual double operator()(const BinnedObservableId& id, DataType d_type = DataType::VALUE) const = 0;

    /**
     * @brief Returns a scalar quantity for a parameter identified by block and LHA id.
     *
     * This overload is intended for proxies already bound to a specific
     * @ref ParameterType.
     *
     * @param block  Block name.
     * @param id     LHA-like identifier inside the block.
     * @param d_type Requested quantity type.
     * @return Requested scalar quantity.
     */
    virtual scalar_t operator()(const std::string& block, const LhaID& id, DataType d_type=DataType::VALUE) const = 0;

    /**
     * @brief Returns the underlying parameter object for a binned observable entry.
     *
     * This is useful when the caller needs access to the full parameter metadata
     * associated with a binned observable, such as uncertainties or binning info.
     *
     * @param id Binned observable identifier.
     * @return Shared pointer to the corresponding @ref Parameter.
     */
    virtual std::shared_ptr<Parameter> get_obs_param(const BinnedObservableId&) const = 0;
};

#endif