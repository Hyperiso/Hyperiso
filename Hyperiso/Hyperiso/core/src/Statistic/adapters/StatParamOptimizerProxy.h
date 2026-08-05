#ifndef STAT_PARAM_OPTIMIZER_PROXY_H
#define STAT_PARAM_OPTIMIZER_PROXY_H

#include "ParamOptimizerAdapter.h"
#include "IStatParamOptimizerProxy.h"

/**
 * @file StatParamOptimizerProxy.h
 * @brief Statistics-layer adapter for optimized batched parameter updates.
 *
 * This class adapts @ref ParamOptimizerAdapter to the statistics module through
 * the @ref IStatParamOptimizerProxy interface.
 *
 * It provides a simple staging API for:
 * - setting scalar values,
 * - setting full parameter objects,
 * - removing parameters,
 * - committing all pending operations in batch,
 * - clearing the current staging area.
 *
 * The proxy keeps a dedicated BSM optimizer and a second optimizer for
 * SM, FLAVOR, DECAY and WILSON parameters.  This preserves the complete
 * parameter scope during likelihood evaluations.
 *
 * This makes the proxy suitable for repeated parameter mutations in statistical
 * tasks such as scans, fitting, profiling, and uncertainty propagation.
 *
 * @see IStatParamOptimizerProxy
 * @see ParamOptimizerAdapter
 */

/**
 * @class StatParamOptimizerProxy
 * @brief Concrete statistics-layer proxy for staged parameter optimization updates.
 *
 * This class is intentionally thin: all logic is delegated to the internal
 * @ref ParamOptimizerAdapter. Its role is to expose a stable, module-local
 * interface to the statistics code.
 */
class StatParamOptimizerProxy : public IStatParamOptimizerProxy {
public:
    /**
     * @brief Constructs the optimizer proxy.
     *
     * Internally initializes the underlying @ref ParamOptimizerAdapter with
     * the parameter types currently relevant to the statistics layer:
     * - @ref ParameterType::SM
     * - @ref ParameterType::FLAVOR
     * - @ref ParameterType::DECAY
     * - @ref ParameterType::WILSON
     *
     * @note The exact list is implementation-defined and may evolve.
     */
    StatParamOptimizerProxy();
    
    /**
     * @copydoc IStatParamOptimizerProxy::set_value
     */
    void set_value(const ParamId& pid, scalar_t v) override;

    /**
     * @copydoc IStatParamOptimizerProxy::set_param
     */
    void set_param(const ParamId& pid, std::shared_ptr<Parameter> p) override;

    /**
     * @copydoc IStatParamOptimizerProxy::remove
     */
    void remove(const ParamId& pid) override;

    /**
     * @copydoc IStatParamOptimizerProxy::commit
     */
    void commit(bool coalesce = true) override;

    /**
     * @copydoc IStatParamOptimizerProxy::clear
     */
    void clear() override;

private:
    ParamOptimizerAdapter& optimizer_for(const ParamId& pid);

    /// Batched updates for the BSM parameter store.
    ParamOptimizerAdapter poa_bsm;

    /// Batched updates for SM, flavor, decay and Wilson parameter stores.
    ParamOptimizerAdapter poa_standard;

};

#endif