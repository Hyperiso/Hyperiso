#ifndef DEPENDENCIES_HELPER_H
#define DEPENDENCIES_HELPER_H

#include "Include.h"

/**
 * @file DependenciesHelper.h
 * @brief Helper for managing parameter–observable dependencies.
 *
 * This helper class centralizes the information about which input parameters
 * (encoded as ParamId) are allowed to affect a given observable or decay.
 *
 * Internally, the dependencies are stored as a static map from DecayId to
 * a set of parameter identifiers. The public interface exposes convenience
 * functions to:
 *   - retrieve the full set of allowed parameters for an observable, and
 *   - test whether a given parameter is allowed for that observable.
 *
 * This is typically used in parameter scans, fits or Monte Carlo sampling to
 * restrict the set of varied parameters to those that are actually relevant
 * to the observable under consideration.
 */
class DependenciesHelper {
public:
    /**
     * @brief Returns the set of allowed parameters for a given observable.
     *
     * This overload takes an Observables enum value, converts it to the
     * corresponding ObservableId using ObservableMapper::to_id(), and then
     * forwards the request to the ObservableId-based overload.
     *
     * @param id High-level observable identifier.
     * @return Set of ParamId objects that are allowed to enter this observable.
     */
    static std::unordered_set<ParamId> get_allowed_parameters(Observables id);

    /**
     * @brief Tests whether a given parameter is allowed for a given observable.
     *
     * This overload takes an Observables enum value, converts it to the
     * corresponding ObservableId and then forwards the check to the
     * ObservableId-based overload.
     *
     * @param id   High-level observable identifier.
     * @param pid  Parameter identifier to be tested.
     * @return true if @p pid is in the allowed parameter set for @p id, false otherwise.
     */

    static bool is_param_allowed(Observables id, ParamId pid);

    /**
     * @brief Returns the set of allowed parameters for a given observable id.
     *
     * The observable id is first mapped to an underlying decay identifier
     * via DecayMapper::get_decay_id(). If the observable is associated with
     * one of the decays listed in the internal dependency table, the
     * corresponding set of ParamId is returned.
     *
     * If the observable does not belong to any declared decay, an error is
     * reported (e.g. via LOG_ERROR in the implementation).
     *
     * @param id Low-level observable identifier.
     * @return Set of ParamId objects allowed for this observable / decay.
     */
    static std::unordered_set<ParamId> get_allowed_parameters(ObservableId id);

    /**
     * @brief Tests whether a given parameter is allowed for a given observable id.
     *
     * This function retrieves the allowed parameter set for @p id using
     * get_allowed_parameters(ObservableId) and checks whether @p pid is
     * contained in that set.
     *
     * @param id   Low-level observable identifier.
     * @param pid  Parameter identifier to be tested.
     * @return true if @p pid is allowed for @p id, false otherwise.
     */
    static bool is_param_allowed(ObservableId id, ParamId pid);
    

private:
    /**
     * @brief Static lookup table of parameter dependencies per decay channel.
     *
     * The map keys are DecayId values (typically obtained from DecayMapper::to_id),
     * and each entry stores the set of parameters that are allowed to enter the
     * computation of observables built from that decay.
     *
     * This table is defined and initialized in the corresponding
     * DependenciesHelper.cpp file.
     */
    static const std::map<DecayId, std::unordered_set<ParamId>> dep_lists;

    /**
     * @brief Static set of common dependencies to all decay channels.
     *
     * This set is defined and initialized in the corresponding
     * DependenciesHelper.cpp file.
     */
    static const std::unordered_set<ParamId> common_deps;

    /**
     * @brief Static set of Wilson coefficients dependencies for all B meson neutral current decays.
     *
     * This set is defined and initialized in the corresponding
     * DependenciesHelper.cpp file.
     */
    static const std::unordered_set<ParamId> B_nc_wilson_deps;

    /**
     * @brief Static set of Wilson coefficients dependencies for all B meson charged current decays.
     *
     * This set is defined and initialized in the corresponding
     * DependenciesHelper.cpp file.
     */
    static const std::unordered_set<ParamId> B_cc_wilson_deps;
};

#endif