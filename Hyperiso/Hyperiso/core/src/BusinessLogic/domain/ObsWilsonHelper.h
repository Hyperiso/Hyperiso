#ifndef OBSWILSONHELPER_H
#define OBSWILSONHELPER_H

#include "Include.h"
#include "Configs.h"
#include "IObsWilsonBuilder.h"
#include "IWilsonFreezer.h"

/**
 * @file ObsWilsonHelper.h
 * @brief Per-manager helper that manages the Wilson group lifecycle for observable computations.
 *
 * The helper decides which Wilson coefficient groups need to be built, which
 * previously active groups should be frozen, and which frozen groups should be
 * unfrozen. Its state is intentionally instance-local: two observable managers
 * with different Wilson builders, scales, or perturbative orders must not share
 * a global cache keyed only by WGroupId.
 */
class ObsWilsonHelper {
public:
    ObsWilsonHelper() = default;

    /**
     * @brief Compatibility constructor used by a few apps/tests to request a reset.
     *
     * Since the helper is no longer backed by static global state, resetting only
     * clears this instance.
     */
    explicit ObsWilsonHelper(bool reset) { if (reset) clear(); }

    /** @brief Clear the per-helper build/freeze state. */
    void clear();

    /**
     * @brief Build or update Wilson groups required by observables.
     *
     * A group is rebuilt when it has never been built by this helper, or when it
     * was previously built with a different matching scale, hadronic scale, or
     * QCD order. This prevents cross-contamination between configurations.
     */
    void build(WilsonBuildConfig config,
               std::shared_ptr<IObsWilsonBuilder>& wil_builder,
               std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer);

    /**
     * @brief Mark a Wilson group as already available in this helper state.
     *
     * Runtime custom Wilson groups are installed through IObsWilsonBuilder and
     * then marked here so DecayParent::enable() does not try to rebuild them via
     * the static Wilson group definitions.
     */
    void mark_built(WGroupId group, const WilsonBuildConfig& config, bool frozen=false);

    /** @brief Mark a group using a minimal single-group build signature. */
    void mark_built(WGroupId group,
                    double matching_scale,
                    double hadronic_scale,
                    QCDOrder order,
                    bool frozen=false);

private:
    struct BuildSignature {
        double matching_scale {};
        double hadronic_scale {};
        QCDOrder order {QCDOrder::NONE};

        bool matches(const WilsonBuildConfig& config) const {
            return matching_scale == config.matching_scale &&
                   hadronic_scale == config.hadronic_scale &&
                   order == config.order;
        }
    };

    struct GroupState {
        bool frozen {false};
        BuildSignature signature {};
    };

    std::unordered_set<WGroupId> get_all_groups(const std::unordered_set<WGroupId>& needed) const;
    std::unordered_set<WGroupId> update_state(const WilsonBuildConfig& config,
                                              std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer);
    static BuildSignature make_signature(const WilsonBuildConfig& config);

    std::unordered_map<WGroupId, GroupState> state;
};

#endif // OBSWILSONHELPER_H
