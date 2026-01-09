#ifndef DECAYPARENT_H
#define DECAYPARENT_H

#include <map>
#include <string>
#include <chrono>
#include <any>
#include <type_traits>

#include "Include.h"
#include "WilsonInterface.h"
#include "ObsWilsonBuilder.h"
#include "ObsWilsonProxy.h"
#include "ObsWilsonHelper.h"
#include "Math.h"
#include "Configs.h"
#include "DefaultConfig.h"
#include "ObservableValue.h"
#include "ObsParameterProxy.h"
#include "ObsPortsConfig.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration;

/**
 * @file DecayParent.h
 * @brief Base class for observable/decay computation modules.
 *
 * This header defines the common interface and lifecycle for "decay" modules
 * (in the broad sense: any observable computation unit).
 *
 * A decay module typically:
 *  - depends on model parameters (SM/BSM/flavor/...) accessed through @ref ObservablePortsConfig
 *  - depends on Wilson coefficients at matching/hadronic scales (built lazily)
 *  - exposes one or multiple observable(s) through @ref compute_observable
 *
 * The class is designed for the business/observable layer:
 *  - it does NOT directly construct low-level providers
 *  - it receives all dependencies via @ref ObservablePortsConfig (ports/adapters)
 *  - it can be enabled/disabled at runtime
 *
 * Wilson coefficients are built "on demand" when @ref enable is called, via
 * @ref ObsWilsonHelper::build. This allows a global Wilson system to be shared,
 * cached, or partially frozen between different decays.
 *
 * @see ObservablePortsConfig
 * @see ObsWilsonHelper
 * @see IObsWilsonBuilder
 * @see IObsWilsonProxy
 */

/**
 * @class DecayParent
 * @brief Abstract base class for a decay/observable module.
 *
 * ### Responsibilities
 *  - Holds the module identifier (@ref DecayId)
 *  - Stores the requested Wilson build configuration (@ref WilsonBuildConfig)
 *    (matching/hadronic scales + QCD order)
 *  - Owns references/pointers to all required "ports":
 *      - Wilson builder/proxy
 *      - parameter proxies
 *      - QCD proxy
 *      - freezer (to freeze/unfreeze Wilson blocks)
 *  - Implements the lifecycle:
 *      - @ref enable builds/activates the module (build Wilson groups, bind proxy, load params)
 *      - @ref disable disables computation without destroying the object
 *
 * ### Lifecycle
 * Typical usage is:
 * @code
 *   MyDecay d(id, muW, muh, QCDOrder::NLO, ports);
 *   d.enable();                 // builds required Wilson groups and loads parameters
 *   auto obs = d.compute_observable(Observables::SOME_OBS);
 * @endcode
 *
 * Calling @ref enable twice is safe (it becomes a no-op after first enable).
 *
 * ### QCD order and MARTY backend
 * - The decay declares a maximum supported order via @ref max_order (default LO).
 * - The constructor stores the requested order, but it is clamped by @ref check_max_order.
 * - If MARTY is enabled (ports.iobs_use_marty), the effective order is limited to LO.
 *
 * @note The actual Wilson construction and caching policy is handled externally via
 *       @ref ObsWilsonHelper and the @ref IWilsonFreezer port.
 */
class DecayParent  {
protected:
    /// Maximum QCD order supported by this decay module (default: LO).
    QCDOrder max_order = QCDOrder::LO;

    /// Reference to the ports/configuration wiring this decay to the framework.
    ObservablePortsConfig& ports;

    /// Wilson builder used to build/add Wilson groups for this decay.
    std::shared_ptr<IObsWilsonBuilder> w_builder;

    /// Wilson proxy used at compute-time to query coefficients (matching/run).
    std::shared_ptr<IObsWilsonProxy> w_proxy;

    /// Port exposing whether the MARTY backend is active.
    std::shared_ptr<IObsCoreAPI<bool>> use_marty;

    /// Parameter proxy for SM-like quantities used by the decay (may be SM/BSM depending on wiring).
    std::shared_ptr<IObsParameterProxy<ParamId, DataType, std::string, LhaID>> p;

    /// QCD proxy (alpha_s, running masses, constants...).
    std::shared_ptr<IObsQCDProxy> iobs_qcdp;

    /// Freezer to freeze/unfreeze Wilson blocks when they are not needed by active decays.
    std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer;

    /// Wilson build configuration used when enabling this decay (scales, order, groups).
    WilsonBuildConfig w_config {};

    /// Whether the decay is enabled (i.e. Wilson groups built and parameters loaded).
    bool enabled {false};

    /// Unique decay identifier.
    DecayId id;

    /**
     * @brief Clamp a requested QCD order to this decay's supported maximum.
     *
     * If @p order is higher than @ref max_order, a warning is emitted and
     * @ref max_order is returned.
     *
     * @param order Requested order.
     * @return Effective order (clamped).
     */
    QCDOrder check_max_order(QCDOrder order) const;

public:
    virtual ~DecayParent() = default;
    
    /**
     * @brief Construct a decay module and bind it to the observable ports.
     *
     * The constructor does NOT build Wilson groups. That happens in @ref enable.
     * It stores the scales and the requested QCD order (clamped to @ref max_order).
     *
     * @param custom_id      Identifier of the decay module.
     * @param matching_scale Matching scale mu_W used for Wilson matching.
     * @param hadronic_scale Hadronic scale mu_h used for running/evolution.
     * @param order          Requested QCD order (will be clamped to @ref max_order).
     * @param ports          Ports configuration providing dependencies.
     */
    DecayParent(DecayId custom_id, double matching_scale, double hadronic_scale, QCDOrder order, ObservablePortsConfig& ports);
    
    /**
     * @brief Bind a Wilson builder to this decay.
     *
     * This is typically used internally by the constructor (binding ports.iobswb),
     * but can be called explicitly if the builder is swapped/replaced.
     *
     * @param wilson_builder Builder used to build Wilson groups needed by this decay.
     */
    void bind_wilson_builder(std::shared_ptr<IObsWilsonBuilder>& wilson_builder);

    /**
     * @brief Enable the decay module (build Wilson groups and load parameters).
     *
     * Steps:
     *  1) Build/update Wilson groups using @ref ObsWilsonHelper::build with @ref w_config
     *  2) Acquire a Wilson proxy from the builder
     *  3) Set a default basis (currently @ref WilsonBasis::B_STANDARD)
     *  4) Call @ref load_params (implemented by derived class)
     *  5) Mark the decay as enabled
     *
     * Calling enable() multiple times is safe: if already enabled it does nothing.
     */
    void enable();

    /**
     * @brief Disable the decay module.
     *
     * This only flips the enabled flag. It does not automatically unbuild or unfreeze
     * Wilson groups; the global policy is handled by @ref ObsWilsonHelper and the freezer port.
     */
    void disable();

     /**
     * @brief Set the QCD order for this decay (one-shot policy).
     *
     * - If MARTY is enabled and @p new_order > LO, it is downgraded to LO.
     * - If the order was never set (w_config.order == NONE), it is set (and clamped).
     * - Otherwise, a warning is emitted and the order is NOT changed.
     *
     * @param new_order Requested order.
     */
    void set_order(QCDOrder new_order);

    /**
     * @brief Load and cache parameters needed by this decay.
     *
     * Derived classes should read all relevant inputs through the proxies
     * (e.g. @ref p, @ref w_proxy, @ref iobs_qcdp) and store them internally.
     *
     * This is called by @ref enable after Wilson has been built and the proxy is ready.
     */
    virtual void load_params() = 0;

    /**
     * @brief Compute an observable given a public observable enum.
     *
     * @param obs Observable enum id.
     * @return One or more results as @ref ObservableValue.
     *
     * @note Some observables are binned; return multiple values or fill @ref ObservableValue::bin.
     */
    virtual std::vector<ObservableValue> compute_observable(Observables obs) = 0;

    /**
     * @brief Compute an observable given an internal observable identifier.
     *
     * @param obs Internal observable id.
     * @return One or more results as @ref ObservableValue.
     */
    virtual std::vector<ObservableValue> compute_observable(ObservableId obs) = 0;

    /**
     * @brief Set a configuration blob for this decay.
     *
     * This generic entry point allows the framework to pass a type-erased config.
     * Most users should prefer inheriting from @ref DecayParentConfigurable<T>
     * to get a typed config setter.
     *
     * @param cfg Type-erased configuration object.
     * @throws std::bad_any_cast if the derived class expects a different type.
     */
    virtual void set_config(std::any cfg) = 0;

    /**
     * @brief Return the decay identifier.
     */
    DecayId get_id() { return id; };
};

/**
 * @class DecayParentConfigurable
 * @brief Typed configuration adapter for @ref DecayParent using std::any.
 *
 * This helper template implements @ref DecayParent::set_config by attempting
 * to extract a configuration object of type @p T from the provided std::any,
 * then forwarding it to a strongly-typed virtual method @ref set_config_spe.
 *
 * This avoids repeating std::any boilerplate in each derived class.
 *
 * Example:
 * @code
 * struct MyConfig { double x; };
 *
 * class MyDecay : public DecayParentConfigurable<MyConfig> {
 * public:
 *   using DecayParentConfigurable::DecayParentConfigurable;
 *   void set_config_spe(MyConfig cfg) override { myx = cfg.x; }
 *   ...
 * };
 * @endcode
 */
template<typename T>
class DecayParentConfigurable : public DecayParent {
public:
    /**
     * @brief Construct a configurable decay module.
     *
     * @param id             Decay id.
     * @param matching_scale Matching scale mu_W.
     * @param hadronic_scale Hadronic scale mu_h.
     * @param order          Requested QCD order.
     * @param ports          Ports configuration.
     */
    DecayParentConfigurable(DecayId id, double matching_scale, double hadronic_scale,
                            QCDOrder order, ObservablePortsConfig& ports)
        : DecayParent(id, matching_scale, hadronic_scale, order, ports)
    {}
    
    /**
     * @brief Type-erased config setter (final).
     *
     * Accepts:
     *  - T
     *  - const T
     *
     * Forwards to @ref set_config_spe(T).
     *
     * @param cfg Type-erased configuration.
     * @throws std::bad_any_cast if cfg does not contain T / const T.
     */
    void set_config(std::any cfg) final override {
        if (auto p = std::any_cast<T>(&cfg)) { set_config_spe(*p); return; }
        if (auto p = std::any_cast<std::add_const_t<T>>(&cfg)) { set_config_spe(*p); return; }
        throw std::bad_any_cast{};
    }

    /**
     * @brief Strongly-typed configuration setter to implement in derived classes.
     *
     * @param config Parsed configuration object.
     */
    virtual void set_config_spe(T config) = 0;
};

/**
 * @brief Specialization for the default decay config type where configuration is ignored.
 *
 * This specialization keeps the same API but implements @ref set_config as a no-op,
 * enabling decays that do not require configuration while still fitting in a generic
 * "configurable" framework.
 */
template<>
struct DecayParentConfigurable<DecayConfig> : DecayParent {
public:
    DecayParentConfigurable(DecayId id, double matching_scale, double hadronic_scale,
                            QCDOrder order, ObservablePortsConfig& ports)
        : DecayParent(id, matching_scale, hadronic_scale, order, ports)
    {}

    /// No-op configuration setter for default config.
    void set_config(std::any) final override {

    }
};

#endif // DECAYPARENT_H