#ifndef HYPERISO_WILSONMANAGER_H
#define HYPERISO_WILSONMANAGER_H

#include <map>
#include <memory>
#include <string>
#include <stdexcept>
#include <set>

#include "Wilson.h"
#include "WilsonGroup.h"
#include "BWilsonGroup.h"
#include "MemoryManager.h"
#include "QCDHelper.h"
#include "Utils.h"
#include "HyperisoMaster.h"
#include "IParamSetter.h"
#include "IWilsonParameters.h"

/**
 * @file WilsonManager.h
 * @brief High-level orchestration of Wilson coefficient groups (matching + running/hadronic).
 *
 * This module provides a “manager” that ties together:
 *  - a set of coefficient groups (@ref CoefficientGroup),
 *  - the dependency engine / block composer (@ref IBlockComposer),
 *  - access to input Wilson blocks (FWCOEF or equivalent) through a parameter proxy,
 *  - the chosen physics model (SM / SUSY / THDM) and backend selection (builtin vs Marty),
 *  - scale management (matching scale mu_W, hadronic scale mu_h).
 *
 * It solves the practical problem of providing consistent triplets of coefficients:
 *  - SM contribution,
 *  - BSM contribution,
 *  - TOTAL = SM + BSM,
 * at the matching scale and then at hadronic scales, for any supported group/basis/order.
 *
 * Internally it uses the DependentBlock/DependentParameter composition mechanisms:
 * - If coefficients must be computed, it composes dependent parameters that call group computations.
 * - If coefficients are provided by input blocks (e.g. FWCOEF), it composes copy/subtraction rules
 *   so the final storage block always has the expected (SM, BSM, TOTAL) triplets.
 *
 * Typical usage:
 * @code
 *   WilsonPortsConfig cfg{ iblock_c, wilson_proxy, use_marty, has_wilson, model_api, scale_setter };
 *   cfg.build_group = ...; // optional factory hook (builder)
 *
 *   auto mgr = CoefficientManager::Builder(
 *       "SM", groups, muW, muH, "NLO", cfg,
 *       {{Model::SM, std::make_shared<WilsonParameterHelper>(iblock_c)}}
 *   );
 *
 *   auto C7tot = mgr->getFullRunCoefficient("B", "C7", "NLO", ContributionType::TOTAL);
 * @endcode
 *
 * @see CoefficientGroup
 * @see IBlockComposer
 * @see IWilsonParameterHelper
 * @see WilsonBlockNames
 */

/**
 * @brief Aggregates “ports” / APIs that connect the manager to the rest of the framework.
 *
 * This struct is deliberately a “wiring” object: it stores shared services and small
 * pieces of policy needed by @ref CoefficientManager.
 *
 * Required services:
 *  - @ref IBlockComposer : registers composed dependent blocks/parameters.
 *  - @ref IParameterProxy : reads existing parameters/blocks (e.g. user-provided FWCOEF).
 *  - @ref ICoreAPI<bool> use_marty : tells which backend is enabled (builtin vs Marty).
 *  - @ref ICoreAPI<bool> has_wilson : whether input Wilson coefficients exist (FWCOEF-like).
 *  - @ref ICoreAPI<Model> model_api : current physics model.
 *  - @ref IParamSetter<ScaleType> scale_setter_api : sets mu_W / mu_h by switching ScaleType.
 *
 * Optional hook:
 *  - build_group: a factory that can build a group for a given (gid, model, backend, contrib, blockname).
 *    This is used notably to build an “SM-only intermediate group” even when the active model is BSM.
 */
struct WilsonPortsConfig {
    /**
     * @brief Constructs the port bundle.
     *
     * @param iblock_c         Block composer used to register dependent blocks/parameters.
     * @param wilson_proxy     Proxy to read/write Wilson parameters by (block, LhaID).
     * @param use_marty        Runtime flag: whether Marty backend is enabled.
     * @param has_wilson       Runtime flag: whether an input FWCOEF block is present.
     * @param model_api        Runtime access to the currently selected physics model.
     * @param scale_setter_api Runtime setter for switching and setting the active scale.
     */
    WilsonPortsConfig(std::shared_ptr<IBlockComposer> iblock_c, std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy, std::shared_ptr<ICoreAPI<bool>> use_marty, std::shared_ptr<ICoreAPI<bool>> has_wilson, std::shared_ptr<ICoreAPI<Model>> model_api, std::shared_ptr<IParamSetter<ScaleType>> scale_setter_api) :
        iblock_c(iblock_c), wilson_proxy(wilson_proxy), 
        use_marty(use_marty), has_wilson(has_wilson),
        model_api(model_api), scale_setter_api(scale_setter_api) {}

    /// Dependency engine used to compose derived blocks/parameters.
    std::shared_ptr<IBlockComposer> iblock_c;

    /// Read-only access to existing blocks/parameters (FWCOEF, matching blocks, etc.).
    std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy;

    /// Backend selector: true => Marty, false => builtin.
    std::shared_ptr<ICoreAPI<bool>> use_marty;

    /// Whether an input Wilson block exists (FWCOEF-style input present).
    std::shared_ptr<ICoreAPI<bool>> has_wilson;

    /// Current model selector (SM/SUSY/THDM...).
    std::shared_ptr<ICoreAPI<Model>> model_api;

    /// Scale setter used to switch and set mu_W / mu_h.
    std::shared_ptr<IParamSetter<ScaleType>> scale_setter_api;

    /**
     * @brief Optional group builder hook.
     *
     * Signature: (group id, model, use_marty, contribution type, storage block name) -> group instance.
     *
     * When provided, the manager uses it to build:
     *  - SM-only intermediate groups (stored in a block name like "<MATCHING>_SM"),
     *  - possibly other specializations depending on the application.
     *
     * If null, the manager falls back to cloning existing groups via @ref CoefficientGroup::get_sm_group().
     */
    std::function<std::shared_ptr<CoefficientGroup>(WGroupId, Model, bool, ContributionType, std::string)> build_group;
};

/**
 * @brief High-level orchestrator for Wilson coefficient groups.
 *
 * The manager owns a registry of groups (by name) and composes all missing blocks/parameters
 * so that:
 *  - matching coefficients are available in a “final matching block” (GroupMapper::... MATCHING),
 *  - hadronic/run coefficients are available in hadronic blocks (ScaleType::HADRONIC, per basis),
 *  - SM/BSM/TOTAL triplets exist consistently for each member coefficient and QCD order.
 *
 * Key responsibilities:
 *
 * 1) Matching initialization (init_group_matching / init_specific_order_group_matching)
 *    - If there is no input Wilson block:
 *        * compute coefficients (group->init(order)) unless "only_total" is requested,
 *        * if model==SM: enforce BSM=0 and TOTAL=SM,
 *        * else (BSM model): ensure an SM intermediate group exists and copy SM -> final,
 *          then compose missing part (TOTAL or BSM) depending on backend.
 *    - If there is input Wilson block (FWCOEF):
 *        * ensure the SM intermediate block exists (SM group),
 *        * compose SM from FWCOEF when available (or from TOT-BSM if both exist),
 *        * compose BSM and TOTAL into the final block with copy/subtraction rules.
 *
 * 2) Hadronic initialization (init_group_hadronic / init_group_hadronic_all_bases)
 *    - Builds hadronic blocks by reading the matching block items and applying
 *      the group running functions per basis and QCD order.
 *
 * 3) Scale management
 *    - set_matching_scale(mu_W) switches ScaleType::MATCHING and sets EW scale,
 *    - set_hadronic_scale(mu_h) switches ScaleType::HADRONIC and sets hadronic scale.
 *
 * Notes on “blocks”:
 *  - Final matching block name is typically: GroupMapper::str(gid, ScaleType::MATCHING).
 *  - SM intermediate block name is: final + "_SM".
 *  - FWCOEF block is "FWCOEF".
 *
 * Lifetime:
 *  - Destructor calls @ref IBlockComposer::remove_all_composed_blocks() to clear composed state.
 */
class CoefficientManager {
private:
    /// Registered coefficient groups, keyed by group name (and possibly also SM intermediate block names).
    std::map<std::string, std::shared_ptr<CoefficientGroup>> coefficientGroups;

    /// Wiring / service access (composer, proxies, APIs, optional builder hook).
    WilsonPortsConfig ports_config;

    /// Logs a detailed message listing available groups when a group name is missing.
    void throw_no_group_error(const std::string& groupName) const;

public:
    CoefficientManager(WilsonPortsConfig ports_config) : ports_config(ports_config) {}
    CoefficientManager(const CoefficientManager&) = delete;
    CoefficientManager operator=(const CoefficientManager&) = delete;

    /**
     * @brief Factory that builds and initializes a manager from pre-created groups.
     *
     * This method:
     *  - ensures Wilson parameter helpers are initialized (if provided),
     *  - registers all provided groups in the manager,
     *  - sets matching/hadronic scales,
     *  - initializes matching + hadronic blocks for each group up to the requested QCD order.
     *
     * @param model Human-readable model string (currently not the canonical model source; see model_api).
     * @param groups Map of groupName -> group instance.
     * @param mu_W Matching scale.
     * @param mu_h Hadronic scale.
     * @param order Max order string ("LO"/"NLO"/"NNLO").
     * @param portconfig Services and hooks needed by the manager.
     * @param wilson_param_helpers Optional per-model helpers that compose WPARAM blocks/matrices.
     */
    static std::shared_ptr<CoefficientManager> Builder(std::string model, std::map<std::string, std::shared_ptr<CoefficientGroup>> groups, double mu_W, double mu_h, std::string order, WilsonPortsConfig portconfig, std::map<Model, std::shared_ptr<IWilsonParameterHelper>> wilson_param_helpers = {});

    /// Switches to matching scale context and sets mu_W via scale_setter_api.
    void set_matching_scale(double mu_W);

    /// Switches to hadronic scale context and sets mu_h via scale_setter_api.
    void set_hadronic_scale(double mu_h);

    /**
     * @brief Ensures an SM intermediate matching block exists and copies SM entries into the final block.
     *
     * This is a crucial step for BSM models: we compute or retrieve SM contributions in a dedicated
     * "<MATCHING>_SM" block, then compose dependent parameters that copy those values into the final block
     * (SM slot of the SM/BSM/TOTAL triplet).
     *
     * @param groupName Group string key.
     * @param order Max order to initialize/copy.
     * @param ports_config Service wiring (uses build_group hook if present).
     */
    void ensure_sm_intermediate_and_copy_to_final(
        const std::string& groupName,
        QCDOrder order,
        WilsonPortsConfig& ports_config
    );

    /**
     * @brief For a given group/order, composes the “missing piece” of the triplet using algebraic relations.
     *
     * Depending on backend and what is computed:
     *  - builtin path tends to compute SM and BSM separately, then compose TOTAL = SM + BSM,
     *  - Marty path tends to provide TOTAL, then we compose BSM = TOTAL - SM (or vice-versa).
     */
    void compose_missing_from_calculation(
        const std::string& groupName,
        QCDOrder order,
        WilsonPortsConfig& ports_config
    );

    /**
     * @brief Composes matching triplets from an input FWCOEF-like block into final storage block.
     *
     * Handles cases where FWCOEF provides:
     *  - SM directly,
     *  - (TOTAL and BSM), allowing reconstruction of SM = TOTAL - BSM,
     *  - only TOTAL or only BSM (fills missing with 0 or with algebra when possible).
     */
    void compose_from_fwcoef(
        const std::string& groupName,
        QCDOrder order,
        WilsonPortsConfig& ports_config
    );

    /// Ensures that, for a given order, SM/BSM/TOTAL slots exist and are set to 0 in the group matching block.
    void ensure_matching_triplet_zeroed(
        const std::string& groupName,
        QCDOrder o
    );

    /// Ensures SM-model semantics: BSM=0 and TOTAL=SM for all coefficients up to max_order in final matching block.
    void ensure_sm_model_triplet_in_matching(
        const std::string& groupName,
        QCDOrder max_order
    );

    /// Ensures final matching block has triplet parameters defined (SM/BSM/TOTAL) defaulting to zero if absent.
    void ensure_final_triplet_defaults_zero(
        const std::string& groupName,
        QCDOrder max_order
    );

    /// Registers a coefficient group under a given name.
    void registerCoefficientGroup(const std::string& groupName, std::shared_ptr<CoefficientGroup> group);

    /// Initializes matching-side parameters for a group up to @p order.
    void init_group_matching(const std::string& groupName, const std::string& order);

    /**
     * @brief Initializes hadronic (running) block for one basis.
     *
     * Reads matching coefficients from the group matching block, applies the basis running function
     * for each contribution type (SM/BSM/TOTAL) and for each order <= requested order, and stores
     * results in the hadronic block (ScaleType::HADRONIC, basis).
     */
    void init_group_hadronic(const std::string& groupName, const std::string& order, WilsonBasis id);

    /// Initializes hadronic blocks for all bases supported by the group.
    void init_group_hadronic_all_bases(const std::string& groupName, const std::string& order);

    /**
     * @brief Initializes matching parameters for a specific QCD order.
     *
     * @param groupName Group string key.
     * @param order QCD order string.
     * @param only_total If true, avoid calling group->init(order) (useful for workflows that only need TOTAL).
     */
    void init_specific_order_group_matching(const std::string& groupName, const std::string& order, bool only_total);

    /// Updates the manager by changing scales (delegates to set_matching_scale / set_hadronic_scale).
    void update(std::string group, double mu_W, double mu_h);
    
    /**
     * @brief Gathers source block names needed to build hadronic/running block dependencies.
     *
     * This collects sources for LO/NLO/NNLO up to the requested order for a given basis.
     */
    void fill_sources_for_group(const std::string & groupName, const std::string& order, std::unordered_map<ParameterType, std::vector<std::string>>& src, WilsonBasis id);

    /**
     * @brief Ensures missing higher-order matching triplets are present/zeroed when only lower order is requested.
     *
     * Example: if order==LO, ensures NLO and NNLO triplets exist and are set to zero, so later computations
     * do not fail when expecting presence of those parameters.
     */
    void fill_matching_groups(const std::string& groupName, const std::string& order);

    /// Returns current model as a string (via model_api).
    std::string getModel();

    /// Retrieves a registered coefficient group by name (throws/logs on missing).
    std::shared_ptr<CoefficientGroup> getCoefficientGroup(const std::string& groupName) const;

    /// Returns all registered groups.
    std::map<std::string, std::shared_ptr<CoefficientGroup>> getGroups();

    /// Returns bases supported by the group (delegates to CoefficientGroup::get_bases()).
    std::unordered_set<WilsonBasis> getGroupBases(WGroupId group);

    /// Returns matching coefficient (raw, order slot) for a given contribution type.
    complex_t getMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type);

    /**
     * @brief Returns “full” matching coefficient including perturbative prefactor powers.
     *
     * Uses fact = alpha_s(mu_W)/(4π) (read from WPARAM_MATCH_SM) and sums:
     *   C_full = C_LO + fact*C_NLO + fact^2*C_NNLO (up to requested order).
     */
    complex_t getFullMatchingCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type);

    /// Returns running/hadronic coefficient (raw, order slot) for given basis and contribution type.
    complex_t getRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type, WilsonBasis basis = WilsonBasis::B_STANDARD);

    /**
     * @brief Returns “full” running coefficient including perturbative prefactor powers.
     *
     * Uses fact = alpha_s(mu_b)/(4π) (read from WPARAM_RUN_SM) and sums:
     *   C_full = C_LO + fact*C_NLO + fact^2*C_NNLO (up to requested order).
     */
    complex_t getFullRunCoefficient(const std::string& groupName, const std::string& coeffName, const std::string& order, ContributionType cont_type, WilsonBasis basis = WilsonBasis::B_STANDARD);

    /// Debug helper (currently stubby in implementation).
    void printGroupCoefficients(const std::string& groupName) const;

    /**
     * @brief Destructor.
     *
     * Clears composed blocks/parameters from the dependency engine to avoid stale
     * state leaking across runs.
     */
    ~CoefficientManager();
};


#endif