#ifndef WILSON_H
#define WILSON_H

#include "Include.h"
#include "Math.h"
#include "Utils.h"
#include "Parameter.h"
#include "IParameterProxy.h"
#include "SourcesView.h"

/**
 * @file Wilson.h
 * @brief Domain objects representing Wilson coefficients and their matching metadata.
 *
 * This header provides:
 * - @ref MatchingInfo: a small aggregate describing how to compute a coefficient at a given QCD order,
 *   and which source parameters it depends on.
 * - @ref WilsonCoefficient: an abstract base class representing one Wilson coefficient, with:
 *   - a name,
 *   - a contribution type (SM/BSM/…),
 *   - a storage block name where values are stored in the Parameters system,
 *   - per-order matching metadata (sources, compute function, LHA id).
 *
 * The class is designed to interact with the framework layers:
 * - It can **reference** values in a Parameters scope through @ref IParameterProxy
 *   (typically a @ref ParameterProxy configured with ParameterType::WILSON).
 * - It can be used with the dependency system: each (order -> MatchingInfo) stores
 *   the set of @ref ParamId sources that should trigger updates.
 *
 * Conventions used:
 * - QCD orders are represented by @ref QCDOrder (LO/NLO/NNLO).
 * - Contribution types are represented by @ref ContributionType (SM/BSM/…).
 * - LHA identifiers are stored as @ref LhaID and follow the mapping conventions of @ref WCoefMapper.
 *
 * @see IParameterProxy
 * @see Parameter
 * @see ParamId
 * @see SourcesView / ParamSrc
 * @see WCoefMapper
 * @see OrderMapper
 */

/**
 * @struct MatchingInfo
 * @ingroup DomainModule
 * @brief Matching metadata for a Wilson coefficient at a specific QCD order.
 *
 * A @ref MatchingInfo describes:
 * - **sources**: the set of parameters required to compute the coefficient at this order,
 * - **compute**: a function taking a @ref ParamSrc view and producing the scalar coefficient value,
 * - **lhaid**: the LHA identifier under which this coefficient/order/type is stored.
 *
 * Default behavior:
 * - `compute` defaults to a lambda returning `0.0` (safe no-op),
 * - `sources` defaults to empty,
 * - `lhaid` is default-constructed unless provided.
 *
 * Notes:
 * - The compute function signature matches the dependency system’s “source view” pattern:
 *   the caller provides a @ref ParamSrc (or similar) that can query source parameters efficiently.
 */
struct MatchingInfo {
    /// Set of source parameter identifiers used by compute().
    std::unordered_set<ParamId> sources = {};

    /**
     * @brief Computation function for this coefficient/order.
     *
     * Takes a @ref ParamSrc providing access to source values and returns
     * the resulting Wilson coefficient (scalar part). Complex coefficients
     * can be represented upstream by multiple scalar entries (e.g. real/imag).
     */
    std::function<scalar_t(const ParamSrc&)> compute =
        [](const auto&) { return 0.0; };

    /// LHA identifier used to store the coefficient value in its storage block.
    LhaID lhaid;

    /// Default ctor: sets compute to a safe no-op (returns 0.0).
    MatchingInfo()
    : compute([](const auto&) { return 0.0; }) {}

    /**
     * @brief Constructs a MatchingInfo with only an LHA id.
     * @param id LHA identifier to associate to this order.
     */
    MatchingInfo(const LhaID& id) : lhaid(id) {}

    /**
     * @brief Full constructor.
     * @param src  Source parameter IDs.
     * @param comp Compute function.
     * @param id   LHA identifier associated with the result.
     */
    MatchingInfo(std::unordered_set<ParamId> src,
                 std::function<scalar_t(const ParamSrc&)> comp,
                 LhaID id)
        : sources(std::move(src)), compute(std::move(comp)), lhaid(std::move(id)) {}

    MatchingInfo(MatchingInfo&& mi) = default;
    MatchingInfo(const MatchingInfo& mi) = default;

    /**
     * @brief Copy assignment operator.
     *
     * Copies sources/compute/lhaid from @p mi.
     */
    MatchingInfo operator=(const MatchingInfo& mi) {
        this->sources = mi.sources;
        this->compute = mi.compute;
        this->lhaid = mi.lhaid;
        return *this;
    }
};

/**
 * @class WilsonCoefficient
 * @ingroup DomainModule
 * @brief Abstract base class representing a Wilson coefficient and its matching information.
 *
 * A Wilson coefficient is represented by:
 * - A **coefficient name** (`coeffName`) used at the domain level.
 * - A **contribution type** (`type`) (SM/BSM/…) controlling naming and storage conventions.
 * - A **storage block** (`storage_block`) indicating where in the Parameters system the
 *   numeric values are stored (typically a WILSON block).
 * - A map `matching_info[order]` storing @ref MatchingInfo per @ref QCDOrder.
 *
 * Construction behavior (as implemented in Wilson.cpp):
 * - If the name ends with `_THDM` or `_SUSY`, the coefficient is marked as BSM.
 * - For each order in {LO, NLO, NNLO}, an entry is created in `matching_info`.
 *
 * Storage / identification:
 * - The LHA id for a given (order,type) is produced via:
 *   @code
 *     WCoefMapper::flha_full(WCoefMapper::enum_elt(get_base_name()), order, typ);
 *   @endcode
 * - get_lhaid_from_name() rebuilds the LHA id using the base coefficient mapping plus
 *   explicit encoding of order and contribution type in the LhaID parts.
 *
 * Reading values:
 * - get_matching_value(...) uses an @ref IParameterProxy (typically WILSON-scoped) to:
 *   - check existence of the entry `(storage_block, id(order, cont_type))`,
 *   - return 0 if missing,
 *   - otherwise read the scalar value and wrap it into complex_t(real=value, imag=0).
 *
 * Ownership flag:
 * - `is_owned` can be used by higher-level code to denote whether this coefficient is
 *   “owned” by a producer (e.g. computed by MARTY / dependent system) vs. externally provided.
 *
 * Extensibility:
 * - The class is abstract: derived types must implement @ref clone().
 */
class WilsonCoefficient {
public:
    WilsonCoefficient() = default;

    /**
     * @brief Constructs a coefficient from a domain name and a storage block.
     *
     * The implementation:
     * - infers ContributionType::BSM if name ends with `_THDM` or `_SUSY`,
     * - initializes per-order matching info using id(order, type).
     *
     * @param name          Domain name of the coefficient (may include suffixes).
     * @param storage_block Block name where coefficient values are stored.
     */
    WilsonCoefficient(const std::string& name, const std::string& storage_block);

    /**
     * @brief Constructs a coefficient from an LHA base name and explicit contribution type.
     *
     * This constructor expects @p name to encode the base coefficient ID (via its parts).
     * It then builds a “full” LHA identifier by appending:
     * - order index (LO=0, NLO=1, NNLO=2),
     * - contribution type index (SM=0, BSM=1, other=2).
     *
     * @param name          LHA identifier representing the base coefficient.
     * @param storage_block Block name where coefficient values are stored.
     * @param ct            Contribution type to associate to this coefficient.
     */
    WilsonCoefficient(const LhaID &name, const std::string& storage_block, ContributionType ct);


    /**
     * @brief Returns the compute function for a given QCD order.
     *
     * @param order QCD order (LO/NLO/NNLO).
     * @return Stored compute functor.
     */
    std::function<scalar_t(const ParamSrc&)> get_func(QCDOrder order);

    /**
     * @brief Returns the source ParamIds required at a given QCD order.
     *
     * @param order QCD order.
     * @return Set of ParamIds used as sources.
     */
    std::unordered_set<ParamId> get_sources(QCDOrder order);

    /**
     * @brief Returns the storage LhaID associated with a given QCD order.
     *
     * @param order QCD order.
     * @return LHA identifier for that order.
     */
    LhaID get_lhaid(QCDOrder order);

    /**
     * @brief Computes the LhaID directly from the coefficient base name and mapping conventions.
     *
     * This is useful when the internal matching_info map does not fully encode
     * how ids are derived from names, or when ids must be reconstructed.
     *
     * @param order QCD order.
     * @return LHA identifier derived from @ref WCoefMapper and the coefficient base name.
     */
    LhaID get_lhaid_from_name(QCDOrder order);

    /**
     * @brief Returns the base name of the coefficient without BSM suffixes.
     *
     * If coeffName ends with `_THDM` or `_SUSY`, the suffix is stripped.
     *
     * @return Base coefficient name.
     */
    std::string get_base_name() const;
    
    /// Sets the domain name of the coefficient.
    void set_name(std::string name) {this->coeffName = name;}

    /**
     * @brief Sets the ownership flag.
     *
     * If already owned and asked to set owned again, this is a no-op.
     *
     * @param owned New ownership value.
     */
    void set_owned(bool owned);

    /// Sets the storage block where the coefficient value is stored.
    void set_storage_block(std::string block_name);

    /// Sets the contribution type (SM/BSM/…).
    void set_contribution_type(ContributionType type);

    /**
     * @brief Reads the coefficient value from the storage backend for a given order/type.
     *
     * This is a *read* operation through an @ref IParameterProxy (commonly WILSON scope).
     * - If the entry does not exist, returns (0,0).
     * - Otherwise returns complex(value, 0).
     *
     * @param order      Order as string (parsed using @ref OrderMapper).
     * @param cont_type  Contribution type (SM/BSM/…).
     * @param wilson_p   Parameter proxy used to access stored values.
     * @return Complex value (imaginary part is 0 with current implementation).
     */
    complex_t get_matching_value(std::string order, ContributionType cont_type, std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_p) const;

    /// Returns the current full name (may include suffixes).
    std::string get_name() const {return this->coeffName;}

    /// Returns the storage block name.
    std::string get_storage_block() const {return this->storage_block;}

    /// Returns the current contribution type.
    ContributionType get_type() {return this->type;}

    /**
     * @brief Returns the LHA id for a given order and contribution type.
     *
     * Delegates to @ref WCoefMapper using the base name mapping.
     *
     * @param order QCD order.
     * @param typ   Contribution type.
     * @return Full LHA identifier.
     */
    LhaID id(QCDOrder order, ContributionType typ) const;

    /// Equality compares name, contribution type, and ownership flag.
    bool operator==(const WilsonCoefficient& other) const;

    /// Inequality convenience.
    bool operator!=(const WilsonCoefficient& other) const { return !(*this == other); }

    virtual ~WilsonCoefficient() = default;

    /**
     * @brief Virtual clone for polymorphic copy.
     *
     * Derived coefficient types must implement this to support copying
     * through base pointers.
     *
     * @return A new heap-allocated copy of the dynamic type.
     */
    virtual std::shared_ptr<WilsonCoefficient> clone() const = 0;

protected:
    /// Domain name of the coefficient (may include suffixes like _THDM/_SUSY).
    std::string coeffName{};

    /// Contribution type for this coefficient (SM by default, inferred in constructors).
    ContributionType type {ContributionType::SM};

    /// Ownership flag (used by higher-level logic to track "produced" coefficients).
    bool is_owned {false};

    /// Block name where this coefficient is stored in the Parameters system.
    std::string storage_block;

    /**
     * @brief Matching metadata indexed by QCD order.
     *
     * Each entry contains:
     * - sources needed to compute this coefficient/order,
     * - compute function,
     * - the LHA id used for storage.
     */
    std::map<QCDOrder, MatchingInfo> matching_info;
};

#endif // WILSON_H