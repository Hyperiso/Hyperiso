#ifndef WILSON_GROUP_H
#define WILSON_GROUP_H

#include <unordered_map>
#include <optional>

#include "Include.h"
#include "Utils.h"
#include "Wilson.h"
#include "IBlockComposer.h"
#include "ICoreAPI.h"
#include "IMartyWilsonProxy.h"
#include "config.hpp"
#include "InterpretedParam.h"
#include "WilsonGroupAdapterConfig.h"

/**
 * @file WilsonGroup.h
 * @brief Domain container grouping related Wilson coefficients with matching and running access.
 *
 * This header defines:
 * - @ref CoefficientGroupSources: per-(basis,order) metadata describing how to build a *running* block:
 *   - which source blocks are required (across ParameterType scopes),
 *   - a function mapping the current coefficient values + source block values into a new coefficient map.
 * - @ref CoefficientGroup: a polymorphic collection of @ref WilsonCoefficient objects representing a
 *   "group" (e.g. a basis-dependent set of coefficients) and providing:
 *   - initialization of matching coefficients as dependent parameters,
 *   - value retrieval at matching scale (mu_W) and hadronic scale (mu_h),
 *   - optional integration points for MARTY-backed generation via adapters.
 *
 * Architectural role:
 * - CoefficientGroup is a *domain-level façade* combining:
 *   - coefficient-level objects (@ref WilsonCoefficient),
 *   - dependency composition (@ref IBlockComposer) to register dependent parameters,
 *   - storage access through @ref IParameterProxy (WILSON scope).
 *
 * Storage conventions:
 * - Matching-scale coefficients are stored in a block identified by @ref block_name
 *   (typically `GroupMapper::str(groupId, ScaleType::MATCHING)`).
 * - Running (hadronic-scale) coefficients are stored in blocks named using:
 *   `GroupMapper::str(groupId, ScaleType::HADRONIC, basis)`.
 *
 * @see WilsonCoefficient
 * @see MatchingInfo
 * @see IBlockComposer
 * @see IParameterProxy
 * @see GroupMapper, OrderMapper, WCoefMapper
 * @see WilsonGroupAdapterConfig
 */

/**
 * @struct CoefficientGroupSources
 * @ingroup DomainModule
 * @brief Source specification and aggregation function for a coefficient group block.
 *
 * This struct is used to describe, for a given Wilson basis and QCD order,
 * how a block of running coefficients should be produced.
 *
 * Members:
 * - @ref sources: describes which *block names* must be available for each ParameterType scope.
 * - @ref func: transforms:
 *     - the current map of coefficient values per order (`ord -> (WCoefId -> value)`),
 *     - the source block view (@ref BlockSrc),
 *   into a new map (WCoefId -> value) for that (basis, order).
 *
 * The default func returns an empty map.
 */
struct CoefficientGroupSources {
    /// Map of required source blocks for each ParameterType scope.
    std::unordered_map<ParameterType, std::vector<std::string>> sources {};

    /**
     * @brief Computes a set of coefficients from existing values and source blocks.
     *
     * Input:
     * - The first argument contains coefficients already available, grouped by QCD order.
     * - The second argument provides read-access to source blocks (BlockSrc).
     *
     * Output:
     * - A flat map from @ref WCoefId to coefficient scalar value (typically for one order/basis).
     */
    std::function<std::unordered_map<WCoefId, scalar_t>(const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>&, const BlockSrc&)> func =
        [](const auto&, const auto&) { return std::unordered_map<WCoefId, scalar_t>(); };
};


/**
 * @class CoefficientGroup
 * @ingroup DomainModule
 * @brief Polymorphic container of WilsonCoefficient objects representing a coefficient group.
 *
 * CoefficientGroup inherits from:
 * `std::map<std::string, std::shared_ptr<WilsonCoefficient>>`
 * mapping coefficient names (keys) to coefficient objects (values).
 *
 * Core responsibilities:
 * - **Claiming coefficients**: @ref claim_coefficients() updates each coefficient so that:
 *   - it is marked owned,
 *   - its storage block matches the group matching block,
 *   - its contribution type matches the group type (SM/BSM/TOTAL...).
 * - **Initialization of matching coefficients**: @ref init(max_order) registers a dependent
 *   parameter for each coefficient and each order up to @p max_order through @ref IBlockComposer.
 * - **Retrieval**:
 *   - @ref get_matching_coefficient() reads the stored coefficient at the matching scale (mu_W).
 *   - @ref get_running_coefficient() reads the stored coefficient at a hadronic scale block (mu_h),
 *     basis-dependent, returning zero if missing.
 *
 * Adapters:
 * - All IO / integration points (dependency, storage, model selection, MARTY) are carried by
 *   @ref WilsonGroupAdapterConfig.
 *
 * Notes on dependent parameter wiring:
 * - init() registers, for each coefficient C and each order,
 *   a wrapper that calls the coefficient’s order-specific compute function:
 *   @code
 *     dep_param->set_expected( coeff.compute(src) );
 *   @endcode
 * - If no compute function is set (null), expected value is set to 0.
 *
 * Copy semantics:
 * - The copy constructor deep-copies coefficients by calling @ref WilsonCoefficient::clone().
 * - Adapter configuration is copied by value (shared_ptr are copied).
 *
 * Extensibility:
 * - This is an abstract base class; derived groups must implement @ref clone().
 */
class CoefficientGroup : public std::map<std::string, std::shared_ptr<WilsonCoefficient>> {
public:
    /**
     * @brief Constructs an empty group with adapter configuration.
     *
     * @param adapters Bundle of adapters/proxies used by the group.
     */
    CoefficientGroup(WilsonGroupAdapterConfig adapters) : adapters(adapters) {};

    /**
     * @brief Deep-copy constructor.
     *
     * Copies:
     * - coefficients via clone(),
     * - group metadata (type, order, id),
     * - sets matching block name from the group id (ScaleType::MATCHING).
     */
    CoefficientGroup(const CoefficientGroup&);

    /// Move construction (default).
    CoefficientGroup(CoefficientGroup&&) = default;

    /**
     * @brief Constructs from an existing map of coefficients and initializes up to NNLO.
     *
     * Inserts all coefficients from @p coeffs, then calls init() with max order NNLO
     * (unless the order system defines NONE as a sentinel).
     *
     * @param coeffs   Map of coefficient names to coefficient instances.
     * @param adapters Adapter config bundle.
     */
    CoefficientGroup(std::map<std::string, std::shared_ptr<WilsonCoefficient>>& coeffs, WilsonGroupAdapterConfig adapters);

    /**
     * @brief Registers dependent parameters for matching coefficients up to a maximum order.
     *
     * For each order in [LO..max_order] and for each coefficient:
     * - builds a wrapper updating the DependentParameter expected value,
     * - registers the dependency through adapters.iblock_c->compose_parameter(...).
     *
     * @param order Maximum QCD order to initialize (LO/NLO/NNLO).
     */
    void init(QCDOrder order);

    /**
     * @brief Initializes sources to build a full running block.
     *
     * This method is intended to set up @ref sources for running blocks, depending on:
     * - the source blocks available,
     * - a chosen Wilson basis,
     * - optional intermediate/interaction controls,
     * - which contribution types are handled.
     *
     * (Implementation-dependent, declared here as part of the group API.)
     */
    void init_full_running_block(const std::unordered_map<ParameterType, std::vector<std::string>> &source_names, WilsonBasis basis, bool inter, std::vector<ContributionType> type);

    /**
     * @brief Returns the SM-only subgroup if the concrete group supports it.
     *
     * Default implementation raises a LogicError and returns nullptr.
     * Concrete derived groups may override.
     */
    virtual std::shared_ptr<CoefficientGroup> get_sm_group() {LOG_ERROR("LogicError", "cannot get_sm_group for virtual group"); return nullptr;}

    /**
     * @brief Reads a coefficient at the matching scale (mu_W).
     *
     * The coefficient is looked up by name in this map.
     * The value is retrieved through the configured WILSON parameter proxy:
     * @code
     *   coeff->get_matching_value(order, cont_type, adapters.wilson_proxy)
     * @endcode
     *
     * @param coeff     Coefficient name (map key).
     * @param order     Order string ("LO", "NLO", "NNLO", ...).
     * @param cont_type Contribution type to request (SM/BSM/TOTAL...).
     * @return Complex value (current implementation usually imag=0 unless stored otherwise).
     * @throws std::out_of_range if coeff is not present in the group.
     */
    complex_t get_matching_coefficient(std::string coeff, std::string order, ContributionType cont_type) const;

    /**
     * @brief Reads a coefficient at the hadronic scale (mu_h) in a given basis.
     *
     * The value is read from the block:
     * `GroupMapper::str(id, ScaleType::HADRONIC, basis)`
     * using the configured WILSON proxy. If the entry does not exist,
     * returns (0,0).
     *
     * @param coeff     Coefficient name (map key).
     * @param order     Order string ("LO"/"NLO"/"NNLO"...).
     * @param cont_type Contribution type to request.
     * @param basis     Wilson basis of the running block (default: B_STANDARD).
     * @return Complex coefficient value; returns (0,0) if missing.
     * @throws std::out_of_range if coeff is not present in the group.
     */
    complex_t get_running_coefficient(std::string coeff, std::string order, ContributionType cont_type, WilsonBasis basis = WilsonBasis::B_STANDARD) const;

    /// Returns the maximum initialized QCD order for this group.
    QCDOrder get_order();

    /**
     * @brief Returns the current matching storage block name.
     *
     * This corresponds to where matching-scale coefficients are stored.
     */
    std::string get_matching_storage_block() const { return block_name; } //TODO to change

    /// Sets the matching storage block name (used by claim_coefficients()).
    void set_matching_storage_block(std::string);

    /// Returns the contribution type associated with this group (SM/BSM/TOTAL...).
    ContributionType get_type() {return this->wilson_type;}

    /**
     * @brief Returns source blocks for a given basis and order.
     * @param ord   QCD order.
     * @param id    Wilson basis.
     * @return Map (ParameterType -> list of block names).
     */
    std::unordered_map<ParameterType, std::vector<std::string>> get_sources(QCDOrder ord, WilsonBasis id) {return this->sources[id][ord].sources;}

    /**
     * @brief Returns the set of bases configured in @ref sources.
     */
    std::unordered_set<WilsonBasis> get_bases() const {
        return get_keys(this->sources);
    }

    /**
     * @brief Returns the computation function for a given basis and order.
     * @param ord QCD order.
     * @param id  Wilson basis.
     */
    std::function<std::unordered_map<WCoefId, scalar_t>(const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>&, const BlockSrc&)> get_func(QCDOrder ord, WilsonBasis id) {return this->sources[id][ord].func;}

    /// Returns the group identifier.
    WGroupId get_group_id() {return id;}

    /// Sets the group identifier (caller must ensure consistency with block naming).
    void set_group_id(WGroupId gid) { id = gid; }

    /// Sets the group contribution type.
    void set_wilson_type(ContributionType ct) { wilson_type = ct; }

    /**
     * @brief Adds running-block sources for a given basis.
     *
     * @param basis Wilson basis.
     * @param m     Map (QCDOrder -> CoefficientGroupSources) for that basis.
     */
    void add_sources(WilsonBasis basis, const std::map<QCDOrder, CoefficientGroupSources>& m) {
        sources[basis] = m;
    }
    
    /**
     * @brief Polymorphic clone.
     *
     * Concrete derived groups must implement this to allow copying through base pointers.
     */
    virtual std::shared_ptr<CoefficientGroup> clone() const = 0;

    virtual ~CoefficientGroup() = default;

protected:
    /**
     * @brief Claims coefficient objects as owned by this group and sets their storage/type.
     *
     * For each coefficient in the group:
     * - sets owned=true,
     * - sets storage block to @ref get_matching_storage_block(),
     * - sets contribution type to @ref wilson_type.
     */
    void claim_coefficients();

    /// Contribution type for this group (default SM).
    ContributionType wilson_type {ContributionType::SM};

    /// Max initialized order (default LO).
    QCDOrder current_order = QCDOrder::LO;

    /// Identifier of this group (used by GroupMapper to build block names).
    WGroupId id;

    /// Block name for matching-scale storage.
    std::string block_name;

    /**
     * @brief Running-block source specification, keyed by (basis -> order).
     */
    std::map<WilsonBasis, std::map<QCDOrder, CoefficientGroupSources>> sources;

    /**
     * @brief Adapter configuration bundle (proxies, composers, core APIs, MARTY hooks).
     */
    WilsonGroupAdapterConfig adapters;
};

/// Streams all coefficients and standard values (LO/NLO/NNLO at matching + hadronic scales).
std::ostream& operator<<(std::ostream& os, const CoefficientGroup& coeffs);

/// Convenience overload for shared_ptr.
std::ostream& operator<<(std::ostream& os, const std::shared_ptr<CoefficientGroup>& coeffs);

#endif // WILSON_GROUP_H