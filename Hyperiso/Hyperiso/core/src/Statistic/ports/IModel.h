#ifndef IMODEL_H
#define IMODEL_H

#include <cstddef>
#include <vector>
#include <map>

#include "Include.h"
#include "ObservableValue.h"

using Vec = std::vector<double>;

/**
 * @file IModel.h
 * @brief Abstract model interface used by the statistical layer.
 *
 * This header defines the model-facing port consumed by the Monte Carlo,
 * likelihood and profiling components. Implementations are responsible for
 * applying fit parameters and nuisance parameters, computing predictions and
 * exposing observable metadata.
 */

/**
 * @class IModel
 * @brief Interface for a model capable of producing observable predictions.
 *
 * The statistical code interacts with concrete model implementations only
 * through this interface. A model receives two parameter maps:
 * - fit parameters @p p, usually the parameters of interest,
 * - nuisance parameters @p eta, usually auxiliary or external inputs.
 *
 * Implementations may cache intermediate state internally, but calls to
 * @ref predict_optimized should return predictions consistent with the
 * provided parameter maps.
 */
class IModel {
public:
    /// Default virtual destructor.
    virtual ~IModel() = default;

    /**
     * @brief Computes model predictions for a given parameter point.
     *
     * @param p   Map of fit-parameter identifiers to values.
     * @param eta Map of nuisance-parameter identifiers to values.
     *
     * @return Predicted observable values grouped by observable identifier.
     */
    virtual std::map<ObservableId, std::vector<ObservableValue>> predict_optimized(const std::map<ParamId, double>& p, const std::map<ParamId, double>& eta) = 0;

    /**
     * @brief Returns the number of currently active binned observables.
     *
     * @return Number of observable bins exposed by the model.
     */
    virtual std::size_t n_observables() const = 0;

    /**
     * @brief Returns the model parameters required by an observable.
     *
     * @param id Observable identifier.
     *
     * @return Set of parameter identifiers on which the observable depends.
     */
    virtual std::unordered_set<ParamId> get_obs_deps(ObservableId id) = 0;

    /**
     * @brief Returns the identifiers of the currently active observable bins.
     *
     * @return Ordered list of binned observable identifiers.
     */
    virtual std::vector<BinnedObservableId> get_obs_ids() = 0;

    /**
     * @brief Forces computation of the currently configured observables.
     *
     * This method is primarily used to initialize or refresh model-side caches
     * before the statistical manager builds likelihood inputs.
     */
    virtual void compute_observables() const = 0;
};

#endif