#ifndef COEFFICIENT_GROUP_BUILDER_H
#define COEFFICIENT_GROUP_BUILDER_H

#include <memory>

#include "GroupDefinition.h"
#include "WilsonGroup.h"
#include "Wilson.h"
#include "WilsonCoefficientRegistry.h"
#include "GenericWilsonGroup.h"

/**
 * @file CoefficientGroupBuilder.h
 * @brief Builder assembling a CoefficientGroup from a GroupDefinition and a CoefficientRegistry.
 *
 * This header defines @ref CoefficientGroupBuilder, a small orchestrator responsible for:
 * - retrieving a @ref GroupDefinition for a requested @ref WGroupId,
 * - creating a concrete @ref CoefficientGroup instance (default: @ref GenericCoefficientGroup),
 * - instantiating each member coefficient via @ref CoefficientRegistry,
 * - wiring coefficient-group "sources" (dependencies) and running optional setup hooks.
 *
 * Key concepts
 * ------------
 * - @ref GroupDefinition describes *what* belongs to a group:
 *   members (WCoef list), dependency sources per basis/order, and optional model setup hooks.
 * - @ref CoefficientRegistry decides *how* to instantiate each coefficient:
 *   it maps (WCoef, Model, Backend) to a concrete @ref WilsonCoefficient subclass factory.
 * - @ref BuildContext carries runtime choices and adapters (proxies, composers, Marty info).
 *
 * Matching block placeholder
 * --------------------------
 * Sources definitions often refer to the matching block name. To keep group definitions
 * reusable, they may use @ref MATCHING_BLOCK_PLACEHOLDER ("$MATCHING_BLOCK$") which is
 * replaced at build time by the actual matching block name computed for the group.
 *
 * @see GroupDefinition
 * @see GroupDefinitions
 * @see CoefficientRegistry
 * @see GenericCoefficientGroup
 * @see BuildContext
 */

/**
 * @class CoefficientGroupBuilder
 * @ingroup DomainModule
 * @brief Builds a ready-to-use @ref CoefficientGroup instance from a @ref BuildContext.
 *
 * The builder is a pure assembly component:
 * - it does not implement coefficient physics,
 * - it does not compute running itself,
 * - it simply instantiates and wires objects according to definitions and registries.
 */
class CoefficientGroupBuilder {
public:
    /**
     * @brief Constructs the builder with a reference to a coefficient registry.
     *
     * @param reg Registry used to create concrete @ref WilsonCoefficient instances.
     */
    explicit CoefficientGroupBuilder(const CoefficientRegistry& reg) : reg_(reg) {}

    /**
     * @brief Builds the coefficient group described by ctx.group_id.
     *
     * Steps performed:
     * 1) Fetch the group definition via @ref GroupDefinitions::get(ctx.group_id).
     * 2) Create a @ref GenericCoefficientGroup initialized with ctx.adapters.
     * 3) Set group id and contribution type (SM/BSM/TOTAL).
     * 4) Determine the matching block name:
     *    - use GroupMapper::str(def.id, ScaleType::MATCHING) by default,
     *    - or ctx.group_name if provided.
     * 5) Copy and patch dependency sources, replacing @ref MATCHING_BLOCK_PLACEHOLDER
     *    with the actual matching block name.
     * 6) Instantiate each coefficient member via @ref CoefficientRegistry::create()
     *    and insert it under its string name (WCoefMapper::str).
     * 7) Run optional setup hooks registered for ctx.model.
     *
     * @param ctx Build context (model/backend/contribution and adapters).
     * @return A fully constructed coefficient group.
     */
    std::shared_ptr<CoefficientGroup> build(const BuildContext& ctx) const;

private:
    /// Registry used to instantiate coefficient objects.
    const CoefficientRegistry& reg_;
};

#endif // COEFFICIENT_GROUP_BUILDER_H