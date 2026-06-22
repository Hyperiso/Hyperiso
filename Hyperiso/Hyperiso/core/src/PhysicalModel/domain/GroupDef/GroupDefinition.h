#ifndef GROUP_DEFINITION_H
#define GROUP_DEFINITION_H

#include <unordered_map>
#include <map>
#include <vector>
#include <functional>
#include <stdexcept>

#include "WilsonGroup.h"
#include "IMartyWilsonPathProxy.h"

/**
 * @file GroupDefinition.h
 * @brief Declarative registry describing Wilson coefficient groups, their sources, and setup hooks.
 *
 * This header introduces a *declarative layer* used to describe Wilson coefficient groups
 * (e.g. B, B', BScalar, charged currents, meson mixing, K, …) without hard-coding construction
 * logic everywhere.
 *
 * A @ref GroupDefinition contains:
 * - a group id (@ref WGroupId),
 * - the list of coefficients belonging to the group (@ref members),
 * - per-basis/per-order source specifications (see @ref CoefficientGroupSources),
 * - optional common and model-specific "setup hooks" executed during group build time (see @ref SetupHook).
 *
 * It also provides the @ref GroupDefinitions namespace which acts as a registry / factory
 * for built-in group definitions and user-registered custom ones.
 *
 * Key concept: MATCHING block placeholder
 * --------------------------------------
 * Many group definitions are written generically and refer to the matching block by using the
 * string placeholder @ref MATCHING_BLOCK_PLACEHOLDER instead of a concrete block name.
 *
 * A builder can then replace it with the actual block:
 * `GroupMapper::str(ctx.group_id, ScaleType::MATCHING)`
 * when instantiating a concrete @ref CoefficientGroup.
 *
 * Built-in vs custom groups
 * ------------------------
 * - Built-in groups are returned by functions like @ref GroupDefinitions::B(), @ref GroupDefinitions::K(), etc.
 * - Custom groups can be registered at runtime via @ref GroupDefinitions::register_custom().
 *
 * The main lookup entry point is:
 * @ref GroupDefinitions::get(WGroupId)
 * which returns the built-in definition if the id is known by @ref GroupMapper, otherwise
 * falls back to the custom registry.
 *
 * @see CoefficientGroup
 * @see CoefficientGroupSources
 * @see WilsonGroupAdapterConfig
 * @see GroupMapper, WCoefMapper
 */

/**
 * @enum Backend
 * @ingroup DomainModule
 * @brief Identifies the computational backend used to obtain coefficients for a group.
 *
 * Typical usage:
 * - Builtin: use internal C++ computations (e.g. `BCoefficientGroup::base_1_LO_calculation`).
 * - Marty:   use MARTY-generated numerical libraries / interfaces when available.
 *
 * This enum is carried in @ref BuildContext and may drive which sources and hooks are applied.
 */
enum class Backend { Builtin, Marty };

/**
 * @struct BuildContext
 * @ingroup DomainModule
 * @brief Build-time context passed to group setup hooks.
 *
 * This context bundles:
 * - adapters and proxies (see @ref WilsonGroupAdapterConfig),
 * - the active physics @ref Model,
 * - the backend choice (built-in vs MARTY),
 * - the active contribution type (SM/BSM/TOTAL),
 * - the group id and an optional human-readable name.
 *
 * Setup hooks can use this context to compose additional dependent blocks/parameters
 * or to override per-model behaviour.
 *
 * @see GroupDefinition::SetupHook
 */
struct BuildContext {
    /// Adapter bundle used by domain objects (storage proxy, composer, core APIs, optional MARTY proxy).
    WilsonGroupAdapterConfig adapters;

    /// Active model used by the core (SM / SUSY / THDM / MARTY).
    Model model   = Model::SM;

    /// Backend selection used to compute or fetch coefficients.
    Backend backend = Backend::Builtin;

     /// Contribution type used for this build (SM/BSM/TOTAL, etc.).
    ContributionType contrib = ContributionType::SM;

    /// Identifier of the group being built.
    WGroupId group_id{};

    /// Optional name for diagnostics / display.
    std::string group_name = "";

    /// Optional path provider used by MARTY-backed Wilson coefficients.
    std::shared_ptr<IMartyWilsonPathProxy> marty_paths;
};

/**
 * @brief Placeholder string used in GroupDefinition sources to refer to the matching block.
 *
 * Many built-in definitions describe their sources using a placeholder instead of a concrete
 * block name. During instantiation, a builder should replace this placeholder with:
 * `GroupMapper::str(ctx.group_id, ScaleType::MATCHING)`.
 *
 * This keeps definitions reusable and avoids duplicating block naming logic.
 */
inline constexpr const char* MATCHING_BLOCK_PLACEHOLDER = "$MATCHING_BLOCK$";

/**
 * @struct GroupDefinition
 * @ingroup DomainModule
 * @brief Declarative description of a Wilson coefficient group.
 *
 * A group definition specifies:
 * - @ref id : a unique group identifier
 * - @ref members : the list of coefficient enums in this group
 * - @ref sources : per-basis and per-order description of required blocks and computation functions
 * - @ref common_setup / @ref setup : optional common and model-specific setup hooks
 *
 * Sources:
 * - `sources[basis][order]` returns a @ref CoefficientGroupSources describing:
 *   - which blocks must exist (across ParameterType scopes),
 *   - which function computes the running coefficients for that (basis, order).
 *
 * Setup hooks:
 * - Hooks can be model-independent (`common_setup`) or stored per model:
 *   `setup[Model::SM]`, `setup[Model::SUSY]`, etc.
 * - Marty builds also inherit SM hooks unless the same hook is already registered
 *   explicitly for Model::MARTY.
 * - Each hook is invoked with a @ref BuildContext and the concrete @ref CoefficientGroup
 *   being built/initialized.
 * - Hooks are typically used to compose auxiliary blocks (e.g. running matrices) needed by
 *   some groups for a given model/backend.
 */
struct GroupDefinition {
    /// Identifier for this group.
    WGroupId id;

    /// Members (coefficient enums) belonging to this group.
    std::vector<WCoef> members;

    /**
     * @brief Running-block sources and computation functions.
     *
     * Keyed by:
     * - WilsonBasis (B_STANDARD, B_TRADITIONAL, ...)
     * - QCDOrder (LO, NLO, NNLO, ...)
     *
     * The source list may include @ref MATCHING_BLOCK_PLACEHOLDER, which should be replaced
     * at build time with the actual matching block name for the group.
     */
    std::unordered_map<WilsonBasis, std::map<QCDOrder, CoefficientGroupSources>> sources;

     /**
     * @brief Setup hook signature.
     *
     * Hooks can:
     * - compose additional dependent blocks via ctx.adapters.iblock_c->compose_block(...),
     * - add/override sources via grp.add_sources(...),
     * - apply per-model or per-backend customizations.
     *
     * @param ctx  Build-time context (model/backend/adapters/contrib/group_id).
     * @param grp  The concrete coefficient group being configured.
     */
    using SetupHook = std::function<void(const BuildContext&, CoefficientGroup&)>;

    /**
     * @brief Setup hooks to run for every model/backend during group construction.
     *
     * Prefer this for auxiliary blocks that are independent of the coefficient
     * backend, e.g. common running matrices.
     */
    std::vector<SetupHook> common_setup;

    /**
     * @brief Model-specific setup hooks to run during group construction.
     *
     * Example:
     * @code
     *   d.setup[Model::SUSY].push_back(&Setup_BScalar_SUSY_Base1_LO);
     * @endcode
     */
    std::unordered_map<Model, std::vector<SetupHook>> setup;

    /**
     * @brief Returns hooks applicable to @p model, including common hooks and
     *        the explicit Marty -> SM compatibility fallback.
     */
    std::vector<SetupHook> hooks_for(Model model) const;
};

/**
 * @namespace GroupDefinitions
 * @ingroup DomainModule
 * @brief Registry/factory of built-in and custom GroupDefinition objects.
 *
 * This namespace provides:
 * - built-in group accessors (B, BPrime, BScalar, charged currents, mixing, K, …),
 * - a unified lookup function @ref get(WGroupId),
 * - a custom registry to allow user-defined group definitions at runtime.
 *
 * Built-in selection:
 * - @ref get(WGroupId) first attempts to map the id through @ref GroupMapper::enum_of().
 * - If recognized, it dispatches to the corresponding built-in accessor.
 * - Otherwise it falls back to the custom registry:
 *   @ref has_custom(), @ref get_custom().
 *
 * Custom registry:
 * - @ref register_custom() stores the definition keyed by its id.
 * - @ref has_custom() checks for existence.
 * - @ref get_custom() returns the stored definition or throws if missing.
 */
namespace GroupDefinitions {

    /// Built-in group definition accessors.
    const GroupDefinition& B();
    const GroupDefinition& BPrime();
    const GroupDefinition& BScalar();
    const GroupDefinition& CC_bc();
    const GroupDefinition& CC_bu();
    const GroupDefinition& CC_cs();
    const GroupDefinition& CC_cd();
    const GroupDefinition& CC_su();
    const GroupDefinition& CC_du();
    const GroupDefinition& MesonMixing();
    const GroupDefinition& K();

    /**
     * @brief Returns the definition for a given group id.
     *
     * If @p g is a known built-in group, returns the corresponding static definition.
     * Otherwise, attempts to retrieve a registered custom definition.
     *
     * @param g Group identifier.
     * @return Reference to the group definition (built-in or custom).
     * @throws std::runtime_error if the id is neither a known built-in nor a registered custom id.
     */
    const GroupDefinition& get(WGroupId g);

    /**
     * @brief Registers a custom group definition at runtime.
     *
     * The definition is stored by id and can later be retrieved via @ref get_custom()
     * or via the generic @ref get() lookup if the id is not recognized as built-in.
     *
     * @param def Custom group definition to register.
     */
    void register_custom(const GroupDefinition& def);

    /**
     * @brief Checks if a custom group definition is registered.
     * @param id Group identifier.
     * @return true if a custom definition exists for @p id.
     */
    bool has_custom(WGroupId id);

    /**
     * @brief Retrieves a previously registered custom group definition.
     * @param id Group identifier.
     * @return Reference to the custom definition.
     * @throws std::runtime_error if no custom definition exists for @p id.
     */
    const GroupDefinition& get_custom(WGroupId id);

}

#endif