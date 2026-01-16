#ifndef OBSWILSONHELPER_H
#define OBSWILSONHELPER_H

#include "Include.h"
#include "Configs.h"
#include "IObsWilsonBuilder.h"
#include "IWilsonFreezer.h"

/**
 * @file ObsWilsonHelper.h
 * @brief High-level helper to manage Wilson group lifecycle for observable computations.
 *
 * @details
 * ObsWilsonHelper is a *stateful orchestration utility* that sits at the boundary
 * between the observable layer and the Wilson infrastructure.
 *
 * Its purpose is to:
 *  - decide which Wilson coefficient groups need to be (re)built,
 *  - avoid rebuilding groups that are already available,
 *  - freeze Wilson blocks that are no longer needed,
 *  - unfreeze previously frozen groups when they are requested again.
 *
 * This logic is essential for observable-driven workflows, where different
 * observables may require different subsets of Wilson coefficient groups,
 * possibly in multiple successive calls.
 *
 * ---
 * ## Conceptual role
 *
 * The helper maintains a **global static state** tracking, for each Wilson group:
 *  - whether it has already been built,
 *  - whether it is currently frozen.
 *
 * When a new observable requests a set of Wilson groups:
 *
 *  1. Groups that are *not requested anymore* are frozen (to avoid recomputation).
 *  2. Groups that are *newly requested* are scheduled for construction.
 *  3. Groups that were frozen but are requested again are unfrozen.
 *
 * Only the minimal required set of Wilson groups is passed to the builder.
 *
 * ---
 * ## Interaction with other components
 *
 * ObsWilsonHelper coordinates three key actors:
 *
 *  - @ref IObsWilsonBuilder  
 *    Responsible for actually building Wilson coefficient groups.
 *
 *  - @ref IWilsonFreezer  
 *    Responsible for freezing/unfreezing Wilson parameter blocks
 *    at both matching and hadronic scales.
 *
 *  - @ref WilsonBuildConfig  
 *    Describes which Wilson groups are required, at which scales and QCD order.
 *
 * ObsWilsonHelper itself does **not** compute Wilson coefficients;
 * it only manages *when* and *whether* they are built or frozen.
 *
 * ---
 * ## Typical usage
 *
 * @code
 * std::shared_ptr<IObsWilsonBuilder> wil_builder = ...;
 * std::shared_ptr<IWilsonFreezer<WGroupId>> freezer = ...;
 *
 * WilsonBuildConfig cfg;
 * cfg.groups = { GroupMapper::to_id(WGroup::B),
 *                GroupMapper::to_id(WGroup::MESON_MIXING) };
 *
 * ObsWilsonHelper::build(cfg, wil_builder, freezer);
 * @endcode
 *
 * Subsequent calls with different group sets will automatically freeze,
 * unfreeze, or rebuild only what is necessary.
 *
 * ---
 * ## State handling
 *
 * The internal state is stored in a static map:
 *
 * @code
 * static std::unordered_map<WGroupId, bool> state;
 * @endcode
 *
 * where:
 *  - `false` → group is active (built and unfrozen),
 *  - `true`  → group is currently frozen.
 *
 * The optional constructor `ObsWilsonHelper(bool reset)` can be used to
 * explicitly reset this state if needed (e.g. for testing or full reinitialization).
 *
 * ---
 * ## Design notes
 *
 * - This class is intentionally *stateless from the user point of view*,
 *   but *stateful internally*.
 * - It avoids expensive recomputation of Wilson coefficients.
 * - It ensures consistency between the observable layer and the Wilson backend.
 *
 * @see IObsWilsonBuilder
 * @see IWilsonFreezer
 * @see WilsonBuildConfig
 * @see ObsWilsonProxy
 */
class ObsWilsonHelper {
public:
    /**
     * @brief Build (or update) Wilson groups required by observables.
     *
     * This method:
     *  - determines which groups must be newly built,
     *  - freezes groups no longer needed,
     *  - unfreezes groups that are requested again,
     *  - delegates the actual construction to the Wilson builder.
     *
     * If no new group needs to be built, the function returns immediately.
     *
     * @param config         Wilson build configuration (groups, scales, order).
     * @param wil_builder    Observable-level Wilson builder.
     * @param iobs_wfreezer  Freezer used to freeze/unfreeze Wilson blocks.
     */
    static void build(WilsonBuildConfig config, std::shared_ptr<IObsWilsonBuilder>& wil_builder, std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer);

    /**
     * @brief Constructor optionally resetting the internal Wilson group state.
     *
     * @param reset If true, clears the internal group state.
     */
    ObsWilsonHelper(bool reset) {reset ? state = {} : state;}

    /// Default constructor (does not reset state).
    ObsWilsonHelper() = default;
private:
    /**
     * @brief Returns the union of already known groups and newly requested ones.
     *
     * @param needed Groups required by the current observable.
     * @return Union of previous and current group sets.
     */
    static std::unordered_set<WGroupId> get_all_groups(const std::unordered_set<WGroupId>& needed);

    /**
     * @brief Updates the internal state and determines which groups must be built.
     *
     * Groups are:
     *  - frozen if no longer needed,
     *  - scheduled for build if newly requested,
     *  - unfrozen if previously frozen but needed again.
     *
     * @param needed        Groups required by the current observable.
     * @param iobs_wfreezer Wilson freezer.
     * @return Set of groups that must be newly built.
     */

    static std::unordered_set<WGroupId> update_state(const std::unordered_set<WGroupId>& needed, std::shared_ptr<IWilsonFreezer<WGroupId>> iobs_wfreezer);

    /// Internal static state tracking frozen/active Wilson groups.
    static inline std::unordered_map<WGroupId, bool> state;
};

#endif // OBSWILSONHELPER_H
