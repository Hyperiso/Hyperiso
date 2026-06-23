#ifndef OBSERVABLE_INTERFACE_H
#define OBSERVABLE_INTERFACE_H

#include <map>
#include <memory>
#include <vector>
#include <string>
#include <exception>
#include <iostream>
#include <cmath>
#include <utility>

#include "Include.h"
#include "Decays.h"
#include "ObsManager.h"
#include "ObsUseMarty.h"
#include "WilsonFreezer.h"
#include "LambdaDecay.h"

/**
 * @file ObservableInterface.h
 * @brief High-level, user-facing entry point to compute flavor observables.
 *
 * @details
 * `ObservableInterface` is the public API intended for end users (scripts, CLI, bindings)
 * to:
 *  - select which observables to compute,
 *  - optionally attach parameter dependences (for scans/likelihood workflows),
 *  - compute theory predictions (possibly binned),
 *  - access experimental central values / uncertainties stored in `FOBS`,
 *  - set / get input parameters in SM/BSM/FLAVOR/etc. parameter spaces.
 *
 * Internally, this interface wires together:
 *  - a Wilson computation pipeline (via `WilsonBuilder` wrapped into `ObsWilsonBuilder`),
 *  - a set of "ports" (proxies) used by the observable layer:
 *      - `ObsParameterProxy` for SM and FLAVOR parameters,
 *      - `ObsQCDProxy` for QCD running (alphas, running masses, constants),
 *      - `ObsUseMarty` for backend capability checks,
 *      - `WilsonFreezer` to freeze/unfreeze Wilson blocks for non-selected groups,
 *  - and an `ObsManager` that owns the registry of decays and selected observables.
 *
 * The design goal is that users only need this class in typical use.
 *
 * @par Typical usage
 * @code
 * ObservableInterface oi;
 *
 * // Add an observable (choose max QCD order for this decay's computation)
 * oi.add_observable(Observables::BR_B_to_Xs_gamma, QCDOrder::NNLO, add_dependencies=true);
 *
 * // Compute prediction (returns possibly multiple bins)
 * auto vals = oi.compute_observable(Observables::BR_B_to_Xs_gamma);
 *
 * // Experimental numbers stored in FOBS (if available)
 * double exp = oi.get_exp_value(Observables::BR_B_to_Xs_gamma);
 * double err = oi.get_exp_uncertainty(Observables::BR_B_to_Xs_gamma, UncertaintyType::COMBINED);
 * @endcode
 *
 * @par Notes on Wilson backend (MARTY)
 * If the selected backend is MARTY, the Wilson side is restricted to LO QCD in this
 * codebase. Higher requested QCD orders for Wilson coefficients are effectively
 * downgraded to LO (see warnings in the Wilson layer).
 *
 * @see ObsManager
 * @see Observable
 * @see DecayParent
 * @see ObservablePortsConfig
 * @see ObsWilsonBuilder
 * @see WilsonFreezer
 */
class ObservableInterface {
private:
    /// Manager that holds registered decays and currently selected observables.
    std::shared_ptr<ObsManager> manager;

    /// Owned ports configuration shared with the manager/decays.
    std::shared_ptr<ObservablePortsConfig> ports;

public:
    /**
     * @brief Construct a ready-to-use observable interface.
     *
     * @details
     * The constructor builds a default "ports" configuration and instantiates an
     * `ObsManager` with default decays enabled.
     *
     * Concretely, it:
     *  - creates a `WilsonBuilder` (initially empty config),
     *  - wraps it as an `ObsWilsonBuilder`,
     *  - creates parameter proxies:
     *      - SM proxy (`ObsParameterProxy()`),
     *      - FLAVOR proxy (`ObsParameterProxy(ParameterType::FLAVOR)`),
     *  - creates the QCD proxy (`ObsQCDProxy`),
     *  - creates the backend capability API (`ObsUseMarty`),
     *  - creates the Wilson freezer (`WilsonFreezer`) bound to the Wilson builder,
     *  - stores them into an `ObservablePortsConfig`,
     *  - finally constructs the `ObsManager` with those ports.
     *
     * @warning The default decays (and their default matching/hadronic scales) are
     * initialized by `ObsManager` using values read from the SM/FLAVOR inputs.
     * If you override parameters, call `reload_params()` or re-enable observables.
     */
    ObservableInterface();

    /**
     * @brief Register a custom decay implementation.
     *
     * @details
     * This allows users to extend the set of computations by providing their own
     * `DecayParent` subclass (typically implementing `load_params()` and
     * `compute_observable()`).
     *
     * The decay will be bound to the current Wilson builder through the manager.
     *
     * @param id  Unique decay identifier.
     * @param ptr Shared pointer to the decay implementation.
     *
     * @see ObsManager::add_custom_decay
     */
    void add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr);

    /**
     * @brief Register a complete custom decay backed by user lambdas.
     *
     * @details
     * This method is the high-level dynamic extension point for observables:
     *  - registers the custom DecayId and all custom ObservableId values,
     *  - installs any lambda-based custom Wilson groups declared in the config,
     *  - creates a @ref LambdaDecay and adds it to the manager,
     *  - optionally selects all observables immediately.
     *
     * The observable lambdas can then use @ref LambdaDecay::W(),
     * @ref LambdaDecay::SM(), @ref LambdaDecay::FLAVOR() and
     * @ref LambdaDecay::QCD() to access the same services as a compiled decay.
     *
     * @param config Runtime/custom decay configuration.
     * @param add_observables If true, all observables declared in @p config are
     *                        added to the manager immediately.
     * @return Reference to `*this` for chaining.
     *
     * @see LambdaDecayConfig
     * @see LambdaObservableConfig
     */
    ObservableInterface& add_lambda_decay(LambdaDecayConfig config, bool add_observables=true);

    /**
     * @brief Add an observable to the manager (enum API).
     *
     * @param obs             Observable enum.
     * @param order           Maximum QCD order requested for this decay.
     * @param add_dependencies If true, automatically adds the full allowed parameter
     *                         dependence set for this observable (via DependenciesHelper).
     * @return Reference to `*this` for chaining.
     *
     * @note The order is stored at the decay level. If multiple observables from the
     * same decay are added, the first effective order selection rules are handled by
     * `DecayParent::set_order()` (see warnings there).
     */
    ObservableInterface& add_observable(Observables obs, QCDOrder order, bool add_dependencies=false);

    /**
     * @brief Add an observable to the manager (internal id API).
     *
     * @param obs             Observable internal id.
     * @param order           Maximum QCD order requested for this decay.
     * @param add_dependencies If true, automatically adds the full allowed parameter
     *                         dependence set for this observable.
     * @return Reference to `*this` for chaining.
     */
    ObservableInterface& add_observable(ObservableId obs, QCDOrder order, bool add_dependencies=false);

    /**
     * @brief Add a binned observable to the manager (internal id API).
     *
     * @param obs             Observable internal id.
     * @param order           Maximum QCD order requested for this decay.
     * @param add_dependencies If true, automatically adds the full allowed parameter
     *                         dependence set for this observable.
     * @return Reference to `*this` for chaining.
     */
    ObservableInterface& add_observable(BinnedObservableId obs, QCDOrder order, bool add_dependencies=false);

    /**
     * @brief Add multiple observables at once (enum map).
     *
     * @param obss            Map (observable -> QCD order).
     * @param add_dependencies If true, add allowed dependences for each observable.
     */
    void add_observables(std::map<Observables, QCDOrder> obss, bool add_dependencies=false);

    /**
     * @brief Add multiple observables at once (internal id map).
     *
     * @param obss            Map (observable id -> QCD order).
     * @param add_dependencies If true, add allowed dependences for each observable.
     */
    void add_observables(std::map<ObservableId, QCDOrder> obss, bool add_dependencies=false);

    /**
     * @brief Add all observables belonging to a decay.
     *
     * @param decay           Decay enum.
     * @param order           Maximum QCD order requested for this decay.
     * @param add_dependencies If true, add allowed dependences for each observable.
     * @param bin q² bin to attach when the decay is binned; ignored otherwise.
     */
    void add_observables(Decays decay, QCDOrder order, bool add_dependencies=false,
                         std::pair<double, double> bin={1.0, 6.0});

    /**
     * @brief Add all observables belonging to a decay.
     *
     * @param decay           Decay dynamic enum.
     * @param order           Maximum QCD order requested for this decay.
     * @param add_dependencies If true, add allowed dependences for each observable.
     * @param bin q² bin to attach when the decay is binned; ignored otherwise.
     */
    void add_observables(DecayId decay, QCDOrder order, bool add_dependencies=false,
                         std::pair<double, double> bin={1.0, 6.0});

    /**
     * @brief Return whether this decay has at least one observable requiring q² bins.
     */
    bool is_decay_binned(Decays decay) const;

    /**
     * @brief Return whether this decay has at least one observable requiring q² bins.
     */
    bool is_decay_binned(DecayId decay) const;

    /**
     * @brief Return whether a specific observable requires q² bins.
     */
    bool is_observable_binned(Observables obs) const;

    /**
     * @brief Return whether a specific observable requires q² bins.
     */
    bool is_observable_binned(ObservableId obs) const;

     /**
     * @brief Manually add a single parameter dependence to an observable (enum API).
     *
     * @details
     * Dependences are used by downstream tooling (scans, profiling, likelihoods) to
     * know which parameters should trigger recomputation.
     *
     * @param obs Observable enum.
     * @param pid Parameter id to add.
     *
     * @note The manager may ignore unsupported dependences based on allow-lists.
     */
    void add_observable_parameter(Observables obs, ParamId pid);

    /**
     * @brief Manually add a single parameter dependence to an observable (id API).
     *
     * @param obs Observable id.
     * @param pid Parameter id to add.
     */
    void add_observable_parameter(ObservableId obs, ParamId pid);

    /**
     * @brief Manually add multiple parameter dependences to an observable (enum API).
     *
     * @param obs  Observable enum.
     * @param pids Set of parameter ids to add.
     */
    void add_observable_parameters(Observables obs, std::unordered_set<ParamId> pids);

     /**
     * @brief Manually add multiple parameter dependences to an observable (id API).
     *
     * @param obs  Observable id.
     * @param pids Set of parameter ids to add.
     */
    void add_observable_parameters(ObservableId obs, std::unordered_set<ParamId> pids);

    /**
     * @brief Compute a theory prediction for an observable (enum API).
     *
     * @details
     * Returns a list of `ObservableValue`:
     *  - for unbinned observables: typically a vector of size 1,
     *  - for binned observables: one entry per bin, with `ObservableValue::bin` set.
     *
     * This call triggers:
     *  - selecting/enabling the corresponding decay,
     *  - ensuring Wilson blocks are built (for needed groups) and frozen otherwise,
     *  - calling the decay’s `compute_observable()`.
     *
     * @param obs Observable enum.
     * @return Vector of predicted values (possibly binned).
     *
     * @see ObservableValue
     * @see ObsManager::evaluate
     */
    std::vector<ObservableValue> compute_observable(Observables obs) const;

    /**
     * @brief Compute a theory prediction for an observable (id API).
     * @param obs Observable id.
     * @return Vector of predicted values (possibly binned).
     */
    std::vector<ObservableValue> compute_observable(ObservableId obs) const;

    /**
     * @brief Compute a theory prediction for an observable (id API).
     * @param obs Observable id with its bin.
     * @return Predicted value.
     */
    ObservableValue compute_observable(BinnedObservableId obs) const;

    /**
     * @brief Remove an observable from the manager (enum API).
     * @param id Observable enum.
     */
    void remove_observable(Observables id);

    /**
     * @brief Remove an observable from the manager (id API).
     * @param id Observable id.
     */
    void remove_observable(ObservableId id);

    /**
     * @brief Remove several observables (enum set).
     * @param ids Set of observable enums.
     */
    void remove_observables(std::unordered_set<Observables> ids);
    
    /**
     * @brief Remove all observables attached to a given decay.
     * @param dec Decay enum.
     */
    void remove_observables(Decays dec);

    /**
     * @brief Get experimental central value for an observable (enum API).
     *
     * @details
     * Reads from the observable parameter space (typically block `"FOBS"`) using the
     * observable proxy stored in `Observable`.
     *
     * @param id Observable enum.
     * @return Experimental central value (if present, provider-defined otherwise).
     */
    scalar_t get_exp_value(Observables id);

    /**
     * @brief Get experimental central value for an observable (id API).
     * @param id Observable id.
     * @return Experimental central value.
     */
    scalar_t get_exp_value(ObservableId id);

    /**
     * @brief Get experimental central value for an observable (id API).
     * @param id Observable id with bin.
     * @return Experimental central value.
     */
    scalar_t get_exp_value(BinnedObservableId id);

    /**
     * @brief Get experimental central value for an observable (id API).
     * @param id Observable id with bin and experience name.
     * @return Experimental central value.
     */
    scalar_t get_exp_value(ExperimentObs id);
    
    /**
     * @brief Get experimental uncertainty for an observable (enum API).
     *
     * @param id     Observable enum.
     * @param u_type Type of uncertainty requested (combined/stat/syst, etc.).
     * @return Experimental uncertainty (provider-defined).
     */
    scalar_t get_exp_uncertainty(Observables id, UncertaintyType u_type=UncertaintyType::COMBINED);

    /**
     * @brief Get experimental uncertainty for an observable (id API).
     * @param id     Observable id.
     * @param u_type Type of uncertainty requested.
     * @return Experimental uncertainty.
     */
    scalar_t get_exp_uncertainty(ObservableId id, UncertaintyType u_type=UncertaintyType::COMBINED);

    /**
     * @brief Get experimental uncertainty for an observable (id API).
     * @param id     Observable id with bins.
     * @param u_type Type of uncertainty requested.
     * @return Experimental uncertainty.
     */
    scalar_t get_exp_uncertainty(BinnedObservableId id, UncertaintyType u_type=UncertaintyType::COMBINED);

    /**
     * @brief Get experimental uncertainty for an observable (id API).
     * @param id     Observable id with bins and experience name.
     * @param u_type Type of uncertainty requested.
     * @return Experimental uncertainty.
     */
    scalar_t get_exp_uncertainty(ExperimentObs id, UncertaintyType u_type=UncertaintyType::COMBINED);

    /**
     * @brief Return the set of currently registered observables.
     * @return Set of observable ids currently stored in the manager.
     */
    std::vector<BinnedObservableId> get_current_observables();

    /**
     * @brief Compute all currently registered observables.
     *
     * @return Map (observable id -> vector of predicted values).
     * @note Each observable may return multiple entries (binned results).
     */
    std::map<ObservableId, std::vector<ObservableValue>> compute_all();

    /**
     * @brief Set an input parameter in a given parameter space.
     *
     * @details
     * This is a convenience wrapper around the global `Parameters` singleton for
     * interactive use.
     *
     * @param block Parameter block name.
     * @param code  LHA-like code inside the block.
     * @param value New value.
     * @param type  ParameterType namespace (SM/BSM/FLAVOR/WILSON/OBSERVABLE/etc.).
     *
     * @warning After setting parameters, cached decay parameters may need refresh.
     * Use `reload_params()` and/or `enable_obs()` depending on your workflow.
     */
    void set_param(const std::string& block, LhaID code, double value, ParameterType type);

    /**
     * @brief Read an input parameter in a given parameter space.
     *
     * @param block Parameter block name.
     * @param code  LHA-like code inside the block.
     * @param type  ParameterType namespace.
     * @return Parameter value.
     */
    double get_param(const std::string& block, LhaID code, ParameterType type);

    /**
     * @brief Retrieve the full allow-list of parameter dependences for an observable (id API).
     *
     * @param id Observable id.
     * @return Set of allowed parameter ids for this observable.
     */
    std::unordered_set<ParamId> get_all_ops_deps(ObservableId id);

    /**
     * @brief Retrieve the full allow-list of parameter dependences for an observable (enum API).
     *
     * @param id Observable enum.
     * @return Set of allowed parameter ids for this observable.
     */
    std::unordered_set<ParamId> get_all_ops_deps(Observables id);

    /**
     * @brief Reload cached parameters for all registered decays.
     *
     * @details
     * Calls `DecayParent::load_params()` on each decay, refreshing internal cached
     * masses/form factors/etc. from the current parameter providers.
     *
     * This is typically used after changing input parameters via `set_param()`.
     */
    void reload_params();

    /**
     * @brief Force enabling/re-enabling currently registered observables.
     *
     * @details
     * This will trigger decay enabling logic in the manager so that Wilson blocks
     * and decay internals are in a consistent state. Useful after large configuration
     * changes or when adding multiple observables.
     */
    void enable_obs();

    /**
     * @brief Set configuration for a given decay.
     *
     * @details
     * The config is passed as `std::any` and interpreted by the decay implementation
     * (see `DecayParentConfigurable<T>`).
     *
     * @param dec    Decay enum.
     * @param config Any-typed configuration payload.
     *
     * @throws std::bad_any_cast if the decay expects a different config type.
     */
    void set_decay_config(Decays dec, std::any config);

    /**
     * @brief Set the thread option for any decay that supports it.
     *
     * Passing 0 delegates to std::thread::hardware_concurrency for supported decays.
     * Unsupported decays keep their default behavior and emit a warning.
     *
     * @param dec Decay enum.
     * @param n_threads Number of threads.
     */
    void set_decay_threads(Decays dec, size_t n_threads);

    /**
     * @brief Set the thread option for the bkstarll decay
     *
     * Because bkstarll is the more time consuming calculation, multithreading is focused on it to improve calculation speed.
     *
     * @param n_threads Number of threads.
     */
    void set_bkstarll_threads(size_t n_threads);

    /**
     * @brief Set the thread option for the bkll decay
     *
     * Because bkll is one of the most time consuming calculation, multithreading is focused on it to improve calculation speed.
     *
     * @param n_threads Number of threads.
     */
    void set_bkll_threads(size_t n_threads);

    /**
     * @brief Set the thread option for the bsphi decay
     *
     * Because bsphi is one of the most time consuming calculation, multithreading is focused on it to improve calculation speed.
     *
     * @param n_threads Number of threads.
     */
    void set_bsphi_threads(size_t n_threads);

    /**
     * @brief Set the thread option for the Lambda_b -> Lambda l+ l- decay.
     *
     * @param n_threads Number of threads.
     */
    void set_lblll_threads(size_t n_threads);

    /**
     * @brief Access the underlying ports configuration (advanced use).
     *
     * @details
     * This exposes the proxies/builder/freezer used by the observable layer and is
     * mainly intended for advanced workflows (custom decays, diagnostics, profiling).
     *
     * @return Reference to the manager-owned ports configuration.
     */
    ObservablePortsConfig& get_ports();
};

#endif