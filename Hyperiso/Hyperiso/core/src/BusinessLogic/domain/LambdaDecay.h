#ifndef LAMBDA_DECAY_H
#define LAMBDA_DECAY_H

#include <any>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "DecayParent.h"
#include "CustomWilsonLambda.h"
#include "mapper_hub.hpp"

class LambdaDecay;

/**
 * @brief User function used to compute a runtime/custom observable.
 *
 * The first argument is the decay execution context. It exposes the same
 * framework services that a hand-written DecayParent subclass would normally
 * use: Wilson proxy, SM/FLAVOR parameter proxies, QCD proxy, bins and the
 * user-provided configuration blob.
 *
 * The second argument is the dynamic observable id currently being evaluated.
 * The function returns one or more ObservableValue entries. Returning several
 * values is the standard way to implement binned observables.
 */
using LambdaObservableComputer =
    std::function<std::vector<ObservableValue>(LambdaDecay&, ObservableId)>;

/**
 * @brief Optional hook called when the lambda decay is enabled/reloaded.
 *
 * Use this for cache preparation, reading parameters once, precomputing common
 * constants, or validating that required Wilson groups are available.
 */
using LambdaDecayLoadHook = std::function<void(LambdaDecay&)>;

/**
 * @brief Optional hook called from LambdaDecay::set_config().
 *
 * This keeps the usual `set_decay_config(..., std::any)` mechanism available
 * for lambda-based decays while still allowing the user to decide how the
 * payload should be interpreted.
 */
using LambdaDecayConfigureHook = std::function<void(LambdaDecay&, const std::any&)>;

/**
 * @struct LambdaObservableConfig
 * @brief Registration and compute function for one runtime/custom observable.
 *
 * A LambdaDecay must own at least one LambdaObservableConfig. Each observable is
 * registered dynamically through ObservableMapper and attached to the parent
 * dynamic DecayId through DecayMapper.
 */
struct LambdaObservableConfig {
    /// Canonical observable name, for example "BR_MY_DECAY".
    std::string canonical;

    /// Optional alternative names accepted by ObservableMapper::id_of().
    std::vector<std::string> aliases {};

    /// Optional FLHA id used for experimental values / binned encodings.
    std::optional<LhaID> flha = std::nullopt;

    /// Optional dependency list used by scans/profilers/likelihood tooling.
    std::unordered_set<ParamId> dependencies {};

    /// User-provided observable computation.
    LambdaObservableComputer compute {};

    LambdaObservableConfig() = default;

    /**
     * @brief Construct a custom observable from a vector-valued compute lambda.
     */
    LambdaObservableConfig(std::string canonical_, LambdaObservableComputer compute_)
        : canonical(std::move(canonical_)), compute(std::move(compute_)) {}

    /**
     * @brief Convenience factory for a scalar, unbinned observable.
     *
     * The user lambda returns only the central value. The factory wraps it into
     * a single ObservableValue whose id is the dynamic id passed at evaluation.
     */
    static LambdaObservableConfig scalar(
        std::string canonical_,
        std::function<double(LambdaDecay&, ObservableId)> compute_
    ) {
        LambdaObservableConfig out;
        out.canonical = std::move(canonical_);
        out.compute = [fn = std::move(compute_)](LambdaDecay& ctx, ObservableId id) {
            return std::vector<ObservableValue>{ObservableValue(id, fn(ctx, id))};
        };
        return out;
    }

    /**
     * @brief Convenience factory for a scalar binned observable.
     *
     * If bins were registered on the decay through ObsManager::add_obs(BinnedObservableId),
     * the lambda is evaluated once per bin and the returned ObservableValue stores
     * the corresponding bin. Without bins, it behaves like a scalar observable and
     * calls the lambda with the conventional bin {0,0}.
     */
    static LambdaObservableConfig binned_scalar(
        std::string canonical_,
        std::function<double(LambdaDecay&, std::pair<double,double>, ObservableId)> compute_
    );
};

/**
 * @struct LambdaDecayConfig
 * @brief Full runtime definition of a custom decay and its observables.
 *
 * This is the observable-side counterpart of @ref CustomWilsonGroupConfig:
 * it lets users create a new decay without writing a new C++ subclass file.
 * The decay can request builtin Wilson groups, custom lambda-based Wilson
 * groups, and one or more lambda-computed observables.
 */
struct LambdaDecayConfig {
    /// Canonical decay name, for example "MY_DECAY".
    std::string canonical;

    /// Optional alternative decay names accepted by DecayMapper::id_of().
    std::vector<std::string> aliases {};

    /// Matching scale used by the decay's Wilson build config.
    double matching_scale = 81.0;

    /// Hadronic scale used by the decay's Wilson build config.
    double hadronic_scale = 4.8;

    /// Requested QCD order for the decay.
    QCDOrder order = QCDOrder::LO;

    /// Maximum QCD order supported by this custom decay.
    QCDOrder max_order = QCDOrder::NNLO;

    /// Builtin or already-registered Wilson groups needed by the decay.
    std::unordered_set<WGroupId> wilson_groups {};

    /**
     * @brief Custom Wilson groups to install before this decay is evaluated.
     *
     * Each group is appended to the shared WilsonBuilder through
     * IObsWilsonBuilder::add_custom_group(), then marked as already built in
     * ObsWilsonHelper so DecayParent::enable() does not try to rebuild it via
     * the static GroupDefinition path.
     */
    std::vector<CustomWilsonGroupConfig> custom_wilson_groups {};

    /// Custom observables owned by this decay. Must not be empty.
    std::vector<LambdaObservableConfig> observables {};

    /**
     * @brief Propagate custom-Wilson matching dependencies to every observable.
     *
     * When a lambda observable is computed from custom Wilson coefficients, the
     * statistical layer ultimately needs the *input* parameters used by those
     * Wilson lambdas.  If this flag is true, ObservableInterface::add_lambda_decay
     * collects the ParamId sources declared in CustomWilsonCoefficientConfig::matching
     * and attaches them to each observable together with the observable-specific
     * dependencies below.
     *
     * Keep this enabled for most dev/plugin use cases. Disable it only when you
     * want to manually control the dependency set of every observable.
     */
    bool propagate_custom_wilson_dependencies {true};

    /// Optional cache/loading hook called by LambdaDecay::load_params().
    LambdaDecayLoadHook load_params {};

    /// Optional type-erased configuration hook called by LambdaDecay::set_config().
    LambdaDecayConfigureHook configure {};
};

/**
 * @class LambdaDecay
 * @brief DecayParent implementation backed by user-provided lambdas.
 *
 * LambdaDecay is intended for dev workflows, plugins and quick phenomenology
 * extensions where adding a full hand-written DecayParent subclass would be too
 * heavy. It still follows the normal framework lifecycle:
 *
 *  - it owns a DecayId,
 *  - it requests Wilson groups through DecayParent::w_config,
 *  - it is enabled/disabled by ObsManager,
 *  - it computes ObservableValue objects for Observable instances.
 *
 * The class does not bypass the mapper layer. A LambdaDecayConfig must first be
 * registered through ObservableInterface::add_lambda_decay(), which creates the
 * DecayId, ObservableId values and DecayMapper links consistently.
 */
class LambdaDecay : public DecayParent {
public:
    /**
     * @brief Construct a lambda-backed decay from a dynamic id and config.
     *
     * @param id Dynamic decay id already registered in DecayMapper.
     * @param config Runtime decay configuration.
     * @param ports Observable ports shared with ObsManager.
     */
    LambdaDecay(DecayId id, LambdaDecayConfig config, ObservablePortsConfig& ports)
        : DecayParent(id, config.matching_scale, config.hadronic_scale, config.order, ports),
          config_(std::move(config))
    {
        this->max_order = config_.max_order;
        this->w_config.groups = config_.wilson_groups;

        for (const auto& custom_group : config_.custom_wilson_groups) {
            this->w_config.groups.insert(custom_group.group);
        }

        for (const auto& obs : config_.observables) {
            if (!obs.compute) {
                LOG_ERROR("ValueError", "Lambda observable", obs.canonical, "has no compute function.");
            }
            computers_.emplace(ObservableMapper::id_of(obs.canonical), obs.compute);
        }
    }

    /**
     * @brief Run the user-defined loading/cache hook, if any.
     */
    void load_params() override {
        if (config_.load_params) {
            config_.load_params(*this);
        }
    }

    /**
     * @brief Compute a builtin observable by converting it to ObservableId.
     */
    std::vector<ObservableValue> compute_observable(Observables obs) override {
        return compute_observable(ObservableMapper::to_id(obs));
    }

    /**
     * @brief Compute a runtime/custom observable by dispatching to its lambda.
     */
    std::vector<ObservableValue> compute_observable(ObservableId obs) override {
        auto it = computers_.find(obs);
        if (it == computers_.end()) {
            LOG_ERROR("KeyError", "Lambda decay", DecayMapper::str(this->id),
                      "does not contain observable", ObservableMapper::str(obs));
        }
        return it->second(*this, obs);
    }

    /**
     * @brief Store and optionally handle a type-erased user config.
     */
    void set_config(std::any cfg) override {
        runtime_config_ = std::move(cfg);
        if (config_.configure) {
            config_.configure(*this, runtime_config_);
        }
    }

    /** @brief Access the immutable LambdaDecayConfig. */
    const LambdaDecayConfig& config() const { return config_; }

    /** @brief Access the last std::any payload passed to set_config(). */
    const std::any& user_config() const { return runtime_config_; }

    /** @brief True if set_config() has received a payload. */
    bool has_user_config() const { return runtime_config_.has_value(); }

    /** @brief Access the Wilson proxy after enable(). */
    IObsWilsonProxy& W() {
        if (!this->w_proxy) {
            LOG_ERROR("LogicError", "Wilson proxy is not available before LambdaDecay is enabled.");
        }
        return *this->w_proxy;
    }

    /** @brief Access the Wilson proxy as a shared pointer after enable(). */
    std::shared_ptr<IObsWilsonProxy> wilson_proxy() {
        if (!this->w_proxy) {
            LOG_ERROR("LogicError", "Wilson proxy is not available before LambdaDecay is enabled.");
        }
        return this->w_proxy;
    }

    /** @brief Access the SM/BSM parameter proxy used by decays. */
    IObsParameterProxy<ParamId, DataType, std::string, LhaID>& SM() { return *this->p; }

    /** @brief Access the FLAVOR parameter proxy from the ports configuration. */
    IObsParameterProxy<ParamId, DataType, std::string, LhaID>& FLAVOR() { return *this->ports.iobspp_flav; }

    /** @brief Access the QCD proxy used by decays. */
    IObsQCDProxy& QCD() { return *this->iobs_qcdp; }

    /** @brief Access all observable ports for advanced use. */
    ObservablePortsConfig& observable_ports() { return this->ports; }

    /** @brief Return the current optional bin list. */
    std::optional<std::vector<std::pair<double,double>>> current_bins() { return this->get_bins(); }

private:
    LambdaDecayConfig config_;
    std::any runtime_config_;
    std::unordered_map<ObservableId, LambdaObservableComputer> computers_;
};

inline LambdaObservableConfig LambdaObservableConfig::binned_scalar(
    std::string canonical_,
    std::function<double(LambdaDecay&, std::pair<double,double>, ObservableId)> compute_
) {
    LambdaObservableConfig out;
    out.canonical = std::move(canonical_);
    out.compute = [fn = std::move(compute_)](LambdaDecay& ctx, ObservableId id) {
        std::vector<ObservableValue> values;
        auto bins = ctx.current_bins();
        if (!bins.has_value() || bins->empty()) {
            values.emplace_back(id, fn(ctx, {0.0, 0.0}, id));
            return values;
        }

        values.reserve(bins->size());
        for (const auto& bin : *bins) {
            values.emplace_back(id, fn(ctx, bin, id), bin);
        }
        return values;
    };
    return out;
}

#endif // LAMBDA_DECAY_H
