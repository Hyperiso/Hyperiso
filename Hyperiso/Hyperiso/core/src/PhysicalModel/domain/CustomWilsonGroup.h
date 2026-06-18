#ifndef CUSTOM_WILSON_GROUP_H
#define CUSTOM_WILSON_GROUP_H

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
#include <utility>

#include "WilsonGroup.h"

/**
 * @file CustomWilsonGroup.h
 * @brief User-assembled CoefficientGroup with explicit coefficients and running rules.
 *
 * This header defines @ref CustomCoefficientGroup, a convenience concrete implementation
 * of @ref CoefficientGroup intended for:
 * - custom groups not present in the built-in @ref GroupDefinitions,
 * - scripting / user extensions,
 * - tests and validation of running/matching pipelines.
 *
 * CustomCoefficientGroup lets the user:
 * - assign an arbitrary group id and user-facing name,
 * - add coefficients manually,
 * - define (basis, order) dependency sources and running function explicitly,
 * - finalize the group by calling @ref finalize(), which delegates to CoefficientGroup::init().
 *
 * @see CoefficientGroup
 * @see CustomWilson
 * @see GroupDefinitions::register_custom
 */

/**
 * @class CustomCoefficientGroup
 * @ingroup DomainModule
 * @brief Concrete CoefficientGroup built manually by the user.
 *
 * This type is especially useful when introducing a new Wilson basis or a new set of
 * coefficients without touching the built-in static group definitions.
 */
class CustomCoefficientGroup : public CoefficientGroup {
public:
    /**
     * @brief Constructs a custom coefficient group.
     *
     * @param adapters Adapters/proxies required by the base CoefficientGroup.
     * @param gid      Group id used in mappers and block naming.
     * @param name     Human-readable group name.
     * @param type     Contribution type (default SM).
     */
    explicit CustomCoefficientGroup(WilsonGroupAdapterConfig adapters, WGroupId gid, const std::string& name,
                                    ContributionType type = ContributionType::SM)
    : CoefficientGroup(adapters) {
        this->id = gid;
        this->name_grp = name;
        this->wilson_type = type;
    }

    /**
     * @brief Polymorphic clone.
     *
     * @return A heap-allocated deep copy of this group.
     */
    std::shared_ptr<CoefficientGroup> clone() const override {
        return std::make_shared<CustomCoefficientGroup>(*this);
    }

    /**
     * @brief Adds a coefficient to the group (fluent API).
     *
     * The coefficient is inserted using its @ref WilsonCoefficient::get_name() as key.
     *
     * @param coef Coefficient to add.
     * @return Reference to *this.
     */
    CustomCoefficientGroup& add_coefficient(const std::shared_ptr<WilsonCoefficient>& coef) {
        this->insert({coef->get_name(), coef});
        this->add_member_id(WCoefMapper::id_of(coef->get_base_name()));
        return *this;
    }

    /**
     * @brief Adds a coefficient with an explicit dynamic id.
     *
     * Use this overload for runtime/custom coefficients where the string stored in
     * the coefficient object may not be convertible to a legacy enum. The explicit
     * @ref WCoefId is recorded as group membership and is later used by
     * @ref CoefficientManager to compose triplets and hadronic blocks.
     *
     * @param id Dynamic coefficient id.
     * @param coef Coefficient implementation.
     * @return Reference to *this.
     */
    CustomCoefficientGroup& add_coefficient(WCoefId id, const std::shared_ptr<WilsonCoefficient>& coef) {
        this->insert({WCoefMapper::str(id), coef});
        this->add_member_id(std::move(id));
        return *this;
    }

    /**
     * @brief Defines dependency sources and running function for a basis and order.
     *
     * This configures the group-level running machinery (hadronic-scale blocks).
     *
     * @param basis        Wilson basis for which this rule applies.
     * @param order        QCD order this rule applies to.
     * @param source_names Block dependencies used by the running function.
     * @param running_func Function producing hadronic coefficients from matching inputs.
     * @return Reference to *this.
     */
    CustomCoefficientGroup& set_basis_order_sources_and_running(
        WilsonBasis basis,
        QCDOrder order,
        std::unordered_map<ParameterType, std::vector<std::string>> source_names,
        std::function<std::unordered_map<WCoefId, scalar_t>(
            const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>&,
            const BlockSrc&
        )> running_func
    ) {
        sources[basis][order].sources = std::move(source_names);
        sources[basis][order].func    = std::move(running_func);
        if (order > current_order) current_order = order;
        return *this;
    }

    /**
     * @brief Default running function that returns the "best available" matching map.
     *
     * Selects NNLO if present, else NLO, else LO, and returns the corresponding map.
     *
     * @param coef_matching Matching coefficient values per order.
     * @param src           Unused (kept to match running function signature).
     * @return Map of coefficients (WCoefId -> value) for the chosen best order.
     */
    static std::unordered_map<WCoefId, scalar_t> identity_running(
        const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& coef_matching,
        const BlockSrc& /*src*/
    ) {
        QCDOrder best = QCDOrder::LO;
        if (coef_matching.count(QCDOrder::NNLO)) best = QCDOrder::NNLO;
        else if (coef_matching.count(QCDOrder::NLO)) best = QCDOrder::NLO;

        return coef_matching.at(best);
    }

    /**
     * @brief Finalizes the group by composing dependent parameters up to @p order.
     *
     * Internally calls @ref CoefficientGroup::init(order), which:
     * - claims coefficients (storage block + contribution type),
     * - composes dependent parameters for matching values via the configured composer.
     *
     * @param order Highest QCD order to initialize (LO/NLO/NNLO).
     * @return Reference to *this.
     */
    CustomCoefficientGroup& finalize(QCDOrder order) {
        this->init(order);
        return *this;
    }
    
    /// Human-readable group name.
    std::string name_grp;
};

#endif // CUSTOM_WILSON_GROUP_H
