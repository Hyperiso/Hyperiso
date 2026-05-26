#ifndef OBSERVABLE_INTERFACE_PROXY_H
#define OBSERVABLE_INTERFACE_PROXY_H

#include <vector>
#include <stdexcept>

#include "ports/IModel.h"
#include "ObservableInterface.h"
#include "StatParamOptimizerProxy.h"
#include "StatParameterProxy.h"
#include "IStatParamOptimizerProxy.h"

/**
 * @file ObservableInterfaceProxy.h
 * @brief Adapter from ObservableInterface to the statistical model interface.
 *
 * The statistical layer expects an @ref IModel. This proxy wraps an
 * @ref ObservableInterface and applies parameter values through an
 * @ref IStatParamOptimizerProxy before requesting observable predictions.
 */

/**
 * @class ObservableInterfaceProxy
 * @brief Concrete @ref IModel implementation backed by ObservableInterface.
 *
 * The proxy translates statistical fit/nuisance parameter maps into updates
 * of the underlying parameter optimizer proxy, commits those updates, then
 * delegates prediction computation to @ref ObservableInterface.
 *
 * @note The constructor taking explicit fit and nuisance lists stores those
 * lists for compatibility with older call sites. The active prediction path
 * requires an optimizer proxy to apply parameter values.
 */
class ObservableInterfaceProxy final : public IModel {
public:
    /**
     * @brief Constructs the proxy with explicit fit and nuisance parameter lists.
     *
     * @param obs       Observable interface to wrap.
     * @param p_specs   Fit-parameter identifiers handled by the statistical layer.
     * @param eta_specs Nuisance-parameter identifiers handled by the statistical layer.
     */
    ObservableInterfaceProxy(
    std::shared_ptr<ObservableInterface> obs,
    std::vector<ParamId> p_specs,
    std::vector<ParamId> eta_specs);
    
    /**
     * @brief Constructs the proxy with a parameter optimizer backend.
     *
     * @param obs  Observable interface to wrap.
     * @param spop Optimizer proxy used to set and commit parameter values.
     */
    ObservableInterfaceProxy(std::shared_ptr<ObservableInterface> obs, std::shared_ptr<IStatParamOptimizerProxy> spop);

    /**
     * @copydoc IModel::n_observables()
     */
    std::size_t n_observables() const override;

    /**
     * @copydoc IModel::get_obs_deps(ObservableId)
     */
    std::unordered_set<ParamId> get_obs_deps(ObservableId id) override;

    /**
     * @copydoc IModel::get_obs_ids()
     */
    std::vector<BinnedObservableId> get_obs_ids() override;

    /**
     * @copydoc IModel::predict_optimized(const std::map<ParamId, double>&, const std::map<ParamId, double>&)
     *
     * @throws std::runtime_error if no optimizer proxy is available.
     */
    std::map<ObservableId, std::vector<ObservableValue>> predict_optimized(
        const std::map<ParamId, double>& p,
        const std::map<ParamId, double>& eta) override;

    /**
     * @copydoc IModel::compute_observables()
     */
    void compute_observables() const;

private:
    std::shared_ptr<ObservableInterface> oi_;             ///< Wrapped observable interface.
    std::shared_ptr<IStatParamOptimizerProxy> spop_;      ///< Parameter setter/committer used before predictions.
    std::vector<ParamId> p_specs_;                        ///< Stored fit-parameter identifiers.
    std::vector<ParamId> eta_specs_;                      ///< Stored nuisance-parameter identifiers.
};

#endif