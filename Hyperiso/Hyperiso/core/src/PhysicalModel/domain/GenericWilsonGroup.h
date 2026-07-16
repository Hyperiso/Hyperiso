#ifndef GENERIC_WILSON_GROUP_H
#define GENERIC_WILSON_GROUP_H

#include "WilsonGroup.h"

/**
 * @file GenericWilsonGroup.h
 * @brief Generic, copyable implementation of a Wilson CoefficientGroup.
 *
 * This header defines @ref GenericCoefficientGroup, a thin concrete class used as the
 * default runtime representation of a @ref CoefficientGroup built from a @ref GroupDefinition.
 *
 * Rationale
 * ---------
 * @ref CoefficientGroup is an abstract type (cloneable) meant to hold a set of named
 * @ref WilsonCoefficient instances, plus metadata (group id, contribution type, sources, etc.).
 * Many groups do not require any custom behavior beyond what @ref CoefficientGroup provides.
 *
 * @ref GenericCoefficientGroup is therefore used as:
 * - the standard implementation returned by builders (e.g. @ref CoefficientGroupBuilder),
 * - a convenient base for copying/cloning groups,
 * - a way to derive an SM-only view of an existing group (see @ref get_sm_group()).
 *
 * @see CoefficientGroup
 * @see GroupDefinition
 * @see CoefficientGroupBuilder
 */

/**
 * @class GenericCoefficientGroup
 * @ingroup DomainModule
 * @brief Default concrete implementation of @ref CoefficientGroup.
 *
 * Inherits all constructors from @ref CoefficientGroup and implements:
 * - @ref clone(): deep copy semantics through the @ref CoefficientGroup copy constructor,
 * - @ref get_sm_group(): helper returning an SM-typed clone (ContributionType::SM).
 *
 * Note:
 * The deep copy behavior relies on @ref CoefficientGroup's copy constructor, which
 * clones each stored @ref WilsonCoefficient via its virtual `clone()` method.
 */
class GenericCoefficientGroup : public CoefficientGroup {
public:
    /// Inherit CoefficientGroup constructors (notably the adapters-based constructor).
    using CoefficientGroup::CoefficientGroup;

    /**
     * @brief Polymorphic clone.
     *
     * @return A heap-allocated deep copy of this group.
     */
    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<GenericCoefficientGroup>(*this);
    }

    /**
     * @brief Returns an SM-only view of this group.
     *
     * Creates a deep copy of the group and forces its contribution type to SM.
     * This is useful when one wants to reuse the same group definition but restrict
     * interpretation of coefficients to SM contributions.
     *
     * @return A deep copy of this group with ContributionType::SM.
     */
    std::shared_ptr<CoefficientGroup> get_sm_group() override {
        auto g = std::make_shared<GenericCoefficientGroup>(*this);
        g->set_wilson_type(ContributionType::SM);
        return g;
    }
};

#endif // GENERIC_WILSON_GROUP_H