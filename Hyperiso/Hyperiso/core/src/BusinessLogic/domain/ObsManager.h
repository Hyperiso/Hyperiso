#ifndef OBSERVABLEMANAGER_H
#define OBSERVABLEMANAGER_H

#include <memory>
#include <map>
#include <unordered_set>
#include <optional>

#include "ParamID.h"
#include "Observable.h"
#include "Decays.h"
#include "IObsParameterProxy.h"
#include "ObsPortsConfig.h"

/**
 * @file ObsManager.h
 * @brief High-level manager for observable evaluation.
 *
 * ObsManager is the main entry point on the **observable side**:
 * it owns a registry of:
 *  - available *decay engines* (@ref DecayParent instances),
 *  - selected *observables* (@ref Observable instances) bound to those decays.
 *
 * It provides convenient methods to:
 *  - add/remove observables to be evaluated,
 *  - evaluate one observable or all registered ones,
 *  - manage observable dependencies (parameter lists),
 *  - configure decays and reload cached parameters,
 *  - register custom decays.
 *
 * ## Design notes
 * - Each observable is backed by a decay. Dynamic @ref ObservableId values are
 *   routed through @ref DecayMapper::get_decay_id, so builtin observables and
 *   custom observable ids registered in the mapper layer follow the same path.
 * - When evaluating a specific observable, ObsManager **selects** the associated decay:
 *   it enables the needed one and disables the others (see @ref select_decay).
 * - Decays lazily build Wilson groups when enabled (through @ref DecayParent::enable),
 *   using the observable ports (Wilson builder, parameter proxies, freezer, etc.).
 *
 * ## Scales defaulting
 * In the default constructor path (when @p init_default_decays is true), matching and
 * hadronic scales are inferred from SM inputs:
 * - mu_W = 2 * m_W
 * - mu_b = m_b(pole) / 2
 *
 * ## Typical usage
 * @code
 * ObservablePortsConfig ports{ ... };
 * ObsManager mgr(ports);
 *
 * // Add one observable with a chosen QCD order
 * mgr.add_obs(Observables::BR_B_to_Xs_gamma, QCDOrder::NLO, add_deps=true);
 *
 * // Evaluate it (returns possibly binned results)
 * auto values = mgr.evaluate(Observables::BR_B_to_Xs_gamma);
 * for (auto& v : values) {
 *     std::cout << v.id.str() << " = " << v.value << "\n";
 * }
 *
 * // Evaluate all registered observables
 * auto all = mgr.evaluate_all();
 * @endcode
 *
 * @see Observable
 * @see DecayParent
 * @see ObservablePortsConfig
 * @see DependenciesHelper
 */
class ObsManager {
public:
    /**
     * @brief Construct an observable manager with the required ports.
     *
     * If @p init_default_decays is true, the manager registers the library default
     * decays (B, K, D, mixing, ...). Their default matching/hadronic scales are derived
     * from SM parameters accessed through @ref ObservablePortsConfig::iobspp_sm.
     *
     * @param obs_port_conf Ports used by decays/observables (Wilson builder, parameter proxies, QCD proxy, freezer, ...).
     * @param init_default_decays Whether to register built-in decays at construction.
     */
    ObsManager(ObservablePortsConfig obs_port_conf, bool init_default_decays=true);

    /**
     * @brief Add an observable (public enum) to the manager.
     *
     * Internally maps to an @ref ObservableId, finds its decay, sets the decay order, and
     * registers an @ref Observable instance bound to that decay.
     *
     * If @p add_deps is true, the observable is populated with its full allowed dependency set
     * from @ref DependenciesHelper.
     *
     * @param id Observable public enum.
     * @param order Requested QCD order for the underlying decay. May be downgraded by the decay.
     * @param add_deps If true, attaches all allowed parameter dependencies to the observable.
     * @return A copy of the manager (builder-like chaining). Note: this returns by value.
     *
     * @see add_obs(ObservableId, QCDOrder, bool)
     */
    ObsManager add_obs(Observables id, QCDOrder order, bool add_deps=false);

    /**
     * @brief Remove an observable (public enum) from the manager.
     *
     * If the observable is not present, a warning may be logged (non-critical path).
     *
     * @param id Observable public enum.
     * @return A copy of the manager (builder-like chaining). Note: this returns by value.
     */
    ObsManager remove_obs(Observables id);

    /**
     * @brief Add an observable (internal ObservableId) to the manager.
     *
     * Steps:
     *  1. Logs the addition.
     *  2. Resolves the parent decay directly from @p id with @ref DecayMapper::get_decay_id.
     *  3. Finds the registered @ref DecayParent instance by @ref DecayId.
     *  4. Sets the decay QCD order (via @ref DecayParent::set_order).
     *  5. Creates and stores an @ref Observable bound to that decay.
     *  6. Optionally attaches dependencies.
     *
     * This overload does not require @p id to be convertible back to the static
     * @ref Observables enum. It therefore supports runtime/custom observables as long
     * as the mapper layer knows their parent decay and that decay has been registered
     * in the manager.
     *
     * @param id Observable internal id.
     * @param order Requested QCD order for the underlying decay. May be downgraded by the decay.
     * @param add_deps If true, attaches all allowed parameter dependencies to the observable.
     * @return A copy of the manager (builder-like chaining). Note: this returns by value.
     */
    ObsManager add_obs(ObservableId id, QCDOrder order, bool add_deps=false);

    /**
     * @brief Add a binned observable to the manager, or add a bin to an existing observable.
     *
     * The parent decay is resolved from @c id.s through @ref DecayMapper::get_decay_id,
     * which keeps the binned path compatible with dynamic/custom observable ids.
     *
     * @param id Observable internal id with bin information.
     * @param order Requested QCD order for the underlying decay. May be downgraded by the decay.
     * @param add_deps If true, attaches all allowed parameter dependencies to the observable.
     * @return A copy of the manager (builder-like chaining). Note: this returns by value.
     */
    ObsManager add_obs(BinnedObservableId id, QCDOrder order, bool add_deps=false);

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
     * @brief Remove an observable (internal ObservableId) from the manager.
     *
     * If the observable is not present, a warning may be logged (non-critical path).
     *
     * @param id Observable internal id.
     * @return A copy of the manager (builder-like chaining). Note: this returns by value.
     */
    ObsManager remove_obs(ObservableId id);

    /**
     * @brief Evaluate a single observable (public enum).
     *
     * This selects (enables) the observable's decay and disables others
     * (see @ref select_decay), then runs @ref Observable::compute.
     *
     * @param id Observable public enum.
     * @return Vector of computed values. May contain multiple entries for binned observables.
     */
    std::vector<ObservableValue> evaluate(Observables id);

    /**
     * @brief Evaluate a single observable (internal ObservableId).
     *
     * This selects (enables) the observable's decay and disables others
     * (see @ref select_decay), then runs @ref Observable::compute.
     *
     * @param id Observable internal id.
     * @return Vector of computed values. May contain multiple entries for binned observables.
     */
    std::vector<ObservableValue> evaluate(ObservableId id);

    /**
     * @brief Evaluate a single observable (internal ObservableId).
     *
     * This selects (enables) the observable's decay and disables others
     * (see @ref select_decay), then runs @ref Observable::compute.
     *
     * @param id Observable internal id with bin information.
     * @return Computed value with bin and id information inside ObservableValue.
     */
    ObservableValue evaluate(BinnedObservableId id);

    /**
     * @brief Evaluate all currently registered observables.
     *
     * Iterates over @ref obss and calls @ref Observable::compute for each.
     * Note: this does **not** call @ref select_decay beforehand, so if multiple
     * observables belong to different decays, each compute() will enable its decay
     * on demand (depending on decay implementation).
     *
     * @return Map from ObservableId to computed values.
     */
    std::map<ObservableId, std::vector<ObservableValue>> evaluate_all();

    /**
     * @brief Register a custom decay implementation.
     *
     * The decay is rebound to the manager Wilson builder (ports), so it can build
     * Wilson groups consistently with the rest of the system.
     *
     * @param id Decay identifier.
     * @param ptr Shared pointer to the decay implementation.
     */
    void add_custom_decay(DecayId id, std::shared_ptr<DecayParent> ptr);

    /**
     * @brief Add a single dependency (parameter) to an observable (public enum).
     *
     * Builtin observables are validated against @ref DependenciesHelper. Runtime
     * observables do not live in that generated table, so their dependencies are
     * accepted as user-declared metadata.
     *
     * @param id Observable public enum.
     * @param param Parameter identifier to add.
     */
    void add_obs_dep(Observables id, ParamId param);

    /**
     * @brief Add a set of dependencies (parameters) to an observable (public enum).
     *
     * Builtin observables are validated against the generated allowed dependency
     * set. Runtime observables accept the provided set directly.
     *
     * @param id Observable public enum.
     * @param params Set of parameters to add.
     */
    void add_obs_deps(Observables id, std::unordered_set<ParamId> params);

    /**
     * @brief Add a single dependency (parameter) to an observable (ObservableId).
     *
     * Builtin observables are validated against @ref DependenciesHelper. Runtime
     * observables accept this dependency directly because the generated helper
     * cannot know plugin/lambda observables in advance.
     *
     * @param id Observable internal id.
     * @param param Parameter identifier to add.
     */
    void add_obs_dep(ObservableId id, ParamId param);

    /**
     * @brief Add a set of dependencies (parameters) to an observable (ObservableId).
     *
     * Builtin observables are validated against the allowed dependency set.
     * Runtime observables accept the provided set directly.
     *
     * @param id Observable internal id.
     * @param params Set of parameters to add.
     */
    void add_obs_deps(ObservableId id, std::unordered_set<ParamId> params);

    /**
     * @brief Get the full allowed dependency set for a given observable.
     *
     * This first returns the dependency set attached to the observable instance.
     * If none is attached, it falls back to the generated DependenciesHelper map
     * for builtin observables. Runtime observables without declared dependencies
     * therefore safely return an empty set instead of throwing.
     *
     * @param id Observable internal id.
     * @return Set of allowed dependencies.
     */
    std::unordered_set<ParamId> get_all_ops_deps(ObservableId id);

    /**
     * @brief Attach all allowed dependencies to an observable (public enum).
     *
     * Equivalent to calling @ref DependenciesHelper::get_allowed_parameters
     * and adding the returned set to the observable.
     *
     * @param id Observable public enum.
     */
    void add_all_obs_deps(Observables id);

    /**
     * @brief Attach all allowed dependencies to an observable (ObservableId).
     *
     * Equivalent to calling @ref DependenciesHelper::get_allowed_parameters
     * and adding the returned set to the observable.
     *
     * @param id Observable internal id.
     */
    void add_all_obs_deps(ObservableId id);

    /**
     * @brief Get the set of observable ids currently registered in the manager.
     * @return Set of ObservableId keys from @ref obss.
     */
    std::vector<BinnedObservableId> get_current_obss();

    /**
     * @brief Retrieve an observable (public enum).
     * @param id Observable public enum.
     * @return Shared pointer to the stored observable.
     */
    std::shared_ptr<Observable> get_obs(Observables id);

    /**
     * @brief Retrieve an observable (ObservableId).
     * @param id Observable internal id.
     * @return Shared pointer to the stored observable.
     */
    std::shared_ptr<Observable> get_obs(ObservableId id);

    /**
     * @brief Access manager ports (non-const reference).
     * @return Reference to the stored ports configuration.
     */
    ObservablePortsConfig& get_ports() {return obs_port_conf;}

    /**
     * @brief Enable the decay needed by an observable and disable all others.
     *
     * This is used by @ref evaluate to avoid keeping multiple decays enabled at once
     * (which can matter if decays build Wilson groups lazily and cache parameters).
     *
     * The decay is selected through @ref DecayMapper::get_decay_id rather than through
     * a static enum conversion, so this also works for registered custom observables.
     *
     * @param id Observable internal id used to determine the corresponding decay.
     */
    void select_decay(ObservableId id);

    /**
     * @brief Reload cached parameters for all registered decays.
     *
     * Calls @ref DecayParent::load_params on each decay instance.
     * Useful after changing global inputs or parameter modes.
     */
    void reload_params();

    /**
     * @brief Enable all observables' decays sequentially.
     *
     * This is a convenience method to (re)initialize decays for all currently registered observables.
     *
     * Each decay is enabled at most once, even when several observables map to it.
     */
    void enable_obs();

    /**
     * @brief Set a decay configuration object.
     *
     * This forwards a type-erased configuration payload to the corresponding decay via
     * @ref DecayParent::set_config. The decay decides how to interpret the payload.
     *
     * @param dec Decay enum.
     * @param config Type-erased configuration object (std::any).
     */
    void set_decay_config(Decays dec, std::any config);

    /**
     * @brief Set the thread option for a decay that supports parallel cache filling.
     *
     * Passing 0 delegates to std::thread::hardware_concurrency for supported decays.
     * Unsupported decays keep their default behavior and emit a warning.
     *
     * @param dec Decay enum.
     * @param n_threads Number of threads.
     */
    ObsManager set_decay_threads(Decays dec, size_t n_threads);

    /**
     * @brief Set the thread option for the B -> K* l+ l- decay.
     *
     * @param n_threads Number of threads.
     */
    ObsManager set_bkstarll_threads(size_t n_threads);

    /**
     * @brief Set the thread option for the B -> K l+ l- decay.
     *
     * @param n_threads Number of threads.
     */
    ObsManager set_bkll_threads(size_t n_threads);

    /**
     * @brief Set the thread option for the Bs -> phi l+ l- decay.
     *
     * @param n_threads Number of threads.
     */
    ObsManager set_bsphi_threads(size_t n_threads);

private:
    /// Registry of decay engines available to the manager.
    std::unordered_map<DecayId, std::shared_ptr<DecayParent>> decays;

    /// Registry of currently selected observables.
    std::map<ObservableId, std::shared_ptr<Observable>> obss;

    /**
     * @brief Ensure a public enum observable is present in the manager.
     *
     * If missing:
     *  - logs an error if @p critical is true,
     *  - logs a warning otherwise.
     *
     * @param id Observable public enum.
     * @param critical Whether missing observable is considered fatal.
     * @return Converted ObservableId.
     */
    ObservableId ensure_present(Observables id, bool critical=true);

    /**
     * @brief Ensure an observable is present in the manager.
     *
     * If missing:
     *  - logs an error if @p critical is true,
     *  - logs a warning otherwise.
     *
     * @param id Observable internal id.
     * @param critical Whether missing observable is considered fatal.
     * @return The same id (for convenience).
     */
    ObservableId ensure_present(ObservableId id, bool critical=true);

    /// Ports used by all decays/observables built by this manager.
    ObservablePortsConfig obs_port_conf;

    /// Currently selected single decay, if the manager is in single-active-decay mode.
    /// A null value means no decay is active yet, or enable_obs() deliberately enabled several decays.
    std::optional<DecayId> active_decay;
};

#endif // OBSERVABLEMANAGER_H
