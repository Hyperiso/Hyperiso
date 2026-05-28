#ifndef HYPERISO_OBSERVABLE_H
#define HYPERISO_OBSERVABLE_H

#include "Include.h"
#include "DecayParent.h"

/**
 * @file Observable.h
 * @brief High-level observable wrapper (user-facing) that connects experimental inputs to theory predictions.
 *
 * An @ref Observable represents one measurable quantity identified by an @ref ObservableId.
 * It provides:
 *  - access to the experimental central value and uncertainties (read from the OBSERVABLE/FOBS block),
 *  - access to the theory prediction computed by its associated @ref DecayParent,
 *  - a dependency list (@ref ParamId) used by higher-level engines (scans, fits, nuisance propagation, caching, etc.).
 *
 * The computation pipeline is:
 *  1. @ref compute() calls @ref DecayParent::enable() (lazy initialization) to ensure required Wilson groups
 *     are built and parameters are loaded (via @ref DecayParent::load_params()).
 *  2. @ref compute() forwards the request to @ref DecayParent::compute_observable(ObservableId).
 *  3. The returned value is a vector of @ref ObservableValue to support:
 *      - single-valued observables,
 *      - binned observables (each entry can carry an optional bin range),
 *      - multi-component outputs when needed.
 *
 * Notes:
 *  - This class is intentionally lightweight: it does not own the physics implementation, only references it.
 *  - Experimental values are fetched through an injected @ref IObsParameterProxy (typically bound to
 *    ParameterType::OBSERVABLE) to keep the observable layer testable and decoupled from the storage backend.
 *
 * Typical usage:
 * @code
 *   // ports.iobspp_obs is usually an ObsParameterProxy(ParameterType::OBSERVABLE)
 *   auto decay = std::make_shared<MyDecay>(DecayId::B_TO_XS_GAMMA, muW, muh, QCDOrder::NNLO, ports);
 *   Observable obs(ObservableId::BR_B_TO_XS_GAMMA, decay, ports.iobspp_obs);
 *
 *   double exp = obs.get_exp_val();
 *   double sig = obs.get_exp_uncertainty(UncertaintyType::COMBINED);
 *   auto theory = obs.compute(); // vector<ObservableValue>
 * @endcode
 *
 * @see DecayParent
 * @see ObservableValue
 * @see IObsParameterProxy
 */
class Observable {
protected:
    /// Identifier of this observable (used for experimental lookup and theory computation routing).
    const ObservableId id;

    /// Physics backend that computes the theory prediction for this observable.
    std::shared_ptr<DecayParent> decay_parent;

    /**
     * @brief Set of parameter identifiers the observable depends on.
     *
     * This is intended for external tooling:
     *  - dependency tracking / invalidation,
     *  - sensitivity studies,
     *  - fit engines (knowing which parameters to vary),
     *  - caching layers.
     */
    std::unordered_set<ParamId> dependences;

    /**
     * @brief Proxy to experimental/observable parameters (typically ParameterType::OBSERVABLE).
     *
     * Expected block usage in current implementation:
     *  - block: "FOBS"
     *  - code : ObservableMapper::flha_of(id)
     */
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_obs;

public:
    /**
     * @brief Construct an Observable.
     *
     * @param id           Observable identifier.
     * @param decay_parent Decay/physics backend used to compute theory prediction.
     * @param iobspp_obs   Proxy used to read experimental values/uncertainties.
     */
    Observable(ObservableId id, std::shared_ptr<DecayParent> decay_parent, std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> iobspp_obs) : id(id), decay_parent(decay_parent), iobspp_obs(iobspp_obs) {}

    /**
     * @brief Return the observable id.
     */
    ObservableId getId() const { return id; }

    /**
     * @brief Get the experimental central value for this observable.
     *
     * Reads from the observable parameter proxy using:
     *  - ParamId(ParameterType::OBSERVABLE, "FOBS", ObservableMapper::flha_of(id))
     *
     * @param bins Experimental bins of the observable
     * @param exp Experiment for the observable
     * @return Experimental central value (as @ref scalar_t).
     */
    scalar_t get_exp_val(std::pair<double, double> bins, std::string exp) const;

    /**
     * @brief Get the experimental uncertainty for this observable.
     *
     * The uncertainty type is mapped to a @ref DataType via @ref UncertaintyTypeMapper::d_type.
     * Typical values include combined/stat/sys depending on your implementation.
     * 
     * @param bins Experimental bins of the observable
     * @param exp Experiment for the observable
     * @param u_type Which uncertainty component to request (default: combined).
     * @return Experimental uncertainty (as @ref scalar_t).
     */
    scalar_t get_exp_uncertainty(std::pair<double, double> bins, std::string exp, UncertaintyType u_type=UncertaintyType::COMBINED) const;

    /**
     * @brief Compute the theory prediction for this observable.
     *
     * This method ensures the underlying @ref DecayParent is enabled before computing.
     * Enabling typically triggers Wilson group building (via ObsWilsonHelper) and parameter loading.
     *
     * @return A vector of observable values (binned or multi-component outputs supported).
     */
    std::vector<ObservableValue> compute() const;

    /**
     * @brief Add a single parameter dependence.
     * @param param_name Parameter id to add to the dependency set.
     */
    void add_dependence(const ParamId& param_name);

    /**
     * @brief Add many parameter dependences.
     *
     * Performs a set union with the existing dependency set.
     *
     * @param param_names Set of parameter ids to add.
     */
    void add_dependences(const std::unordered_set<ParamId>& param_names);

    /**
     * @brief Access the dependency set (read-only).
     * @return Reference to internal dependency set.
     */
    const std::unordered_set<ParamId>& get_dependences() const;
}; 


#endif // HYPERISO_OBSERVABLE_H
