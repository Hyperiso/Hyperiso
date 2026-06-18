#ifndef CUSTOM_WILSON_LAMBDA_H
#define CUSTOM_WILSON_LAMBDA_H

#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "CustomWilson.h"
#include "CustomWilsonGroup.h"
#include "WilsonBlockNames.h"

/**
 * @file CustomWilsonLambda.h
 * @brief High-level configuration helpers for user-defined Wilson groups and coefficients.
 *
 * The lower domain layer already contains @ref CustomWilson and
 * @ref CustomCoefficientGroup. This header adds small value objects that make the
 * public API easy to use from examples, bindings and configuration-driven code:
 *
 *   - @ref CustomWilsonCoefficientConfig describes one coefficient and the lambdas
 *     used to compute its matching values at LO/NLO/NNLO.
 *   - @ref CustomWilsonGroupConfig describes one group, its coefficients, its
 *     matching/running scales and its optional running lambdas.
 *   - @ref make_custom_wilson_group turns the config into a fully wired
 *     @ref CustomCoefficientGroup ready to be registered in @ref CoefficientManager.
 *
 * All identifiers are dynamic ids (@ref WGroupId and @ref WCoefId). Builtin enums
 * are not required, so this path works for symbols registered at runtime through
 * @ref GroupMapper and @ref WCoefMapper.
 */

/**
 * @struct CustomWilsonCoefficientConfig
 * @brief Configuration of one user-defined Wilson coefficient.
 *
 * A coefficient is identified by a dynamic @ref WCoefId and receives one matching
 * rule per QCD order. Each rule is a standard @ref MatchingInfo containing:
 *   - parameter dependencies, as @ref ParamId values,
 *   - a compute lambda of type `scalar_t(const ParamSrc&)`,
 *   - the full FLHA id under which the result is stored.
 *
 * If an order is not configured, a zero-valued default rule is installed for that
 * order. This keeps the downstream matching/hadronic composition robust when a
 * caller requests NLO/NNLO while only LO was supplied.
 */
struct CustomWilsonCoefficientConfig {
    /// Dynamic coefficient id. It must be known by WCoefMapper and have an FLHA base key.
    WCoefId id;

    /// Per-order matching rules.
    std::map<QCDOrder, MatchingInfo> matching;

    CustomWilsonCoefficientConfig() = default;
    explicit CustomWilsonCoefficientConfig(WCoefId id) : id(std::move(id)) {}

    /**
     * @brief Attach a matching lambda to a QCD order.
     *
     * @param order QCD order to configure.
     * @param sources Parameter dependencies needed by @p compute.
     * @param compute User function computing the coefficient from a @ref ParamSrc.
     * @param contribution Contribution slot this rule fills; used only to build the default FLHA id.
     * @return Reference to this config for chaining.
     */
    CustomWilsonCoefficientConfig& set_matching(
        QCDOrder order,
        std::unordered_set<ParamId> sources,
        std::function<scalar_t(const ParamSrc&)> compute,
        ContributionType contribution = ContributionType::SM
    ) {
        matching[order] = MatchingInfo(
            std::move(sources),
            std::move(compute),
            WCoefMapper::flha_full(id, order, contribution)
        );
        return *this;
    }
};

/**
 * @struct CustomWilsonGroupConfig
 * @brief Configuration of a user-defined Wilson group built from lambdas.
 *
 * The group config contains:
 *   - a dynamic group id,
 *   - group-level scales and maximum QCD order,
 *   - the contribution slot computed by coefficient lambdas,
 *   - a list of coefficient configs,
 *   - optional running rules per Wilson basis and QCD order.
 *
 * If no running rule is supplied, an identity running rule is installed for
 * @ref WilsonBasis::B_STANDARD at LO/NLO/NNLO. It forwards the best available
 * matching coefficients to the hadronic block without additional evolution.
 */
struct CustomWilsonGroupConfig {
    /// Dynamic group id. Register it with GroupMapper before use if it is custom.
    WGroupId group;

    /// Optional display name used by CustomCoefficientGroup; the mapper name remains authoritative.
    std::string display_name;

    /// Matching scale used when the group is added to the manager.
    double matching_scale {81.0};

    /// Hadronic/running scale used when the group is added to the manager.
    double hadronic_scale {4.8};

    /// Highest QCD order to initialize.
    QCDOrder order {QCDOrder::LO};

    /// Contribution slot computed by the custom matching lambdas.
    ContributionType contribution {ContributionType::SM};

    /// Coefficients belonging to the group.
    std::vector<CustomWilsonCoefficientConfig> coefficients;

    /// Running rules keyed by basis and QCD order.
    std::unordered_map<WilsonBasis, std::map<QCDOrder, CoefficientGroupSources>> running;

    /// If true, install identity running for B_STANDARD when no running rules were provided.
    bool install_identity_running_if_empty {true};

    CustomWilsonGroupConfig() = default;
    explicit CustomWilsonGroupConfig(WGroupId group) : group(std::move(group)) {}

    /**
     * @brief Add a coefficient config to the group.
     */
    CustomWilsonGroupConfig& add_coefficient(CustomWilsonCoefficientConfig coef) {
        coefficients.emplace_back(std::move(coef));
        return *this;
    }

    /**
     * @brief Attach a running lambda for a basis/order pair.
     *
     * The lambda receives all matching coefficients already available by QCD order
     * and a @ref BlockSrc for additional source blocks. It must return the flat map
     * of hadronic coefficients for the requested order.
     */
    CustomWilsonGroupConfig& set_running(
        WilsonBasis basis,
        QCDOrder order,
        std::unordered_map<ParameterType, std::vector<std::string>> sources,
        std::function<std::unordered_map<WCoefId, scalar_t>(
            const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>&,
            const BlockSrc&
        )> compute
    ) {
        CoefficientGroupSources s;
        s.sources = std::move(sources);
        s.func = std::move(compute);
        running[basis][order] = std::move(s);
        return *this;
    }
};

/**
 * @brief Identity running rule used as a safe default for custom groups.
 *
 * It selects the best matching order available in the input map: NNLO, else NLO,
 * else LO. If no matching coefficient exists, it returns an empty map.
 */
inline std::unordered_map<WCoefId, scalar_t> custom_identity_running(
    const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
    const BlockSrc&
) {
    if (auto it = matching.find(QCDOrder::NNLO); it != matching.end()) return it->second;
    if (auto it = matching.find(QCDOrder::NLO);  it != matching.end()) return it->second;
    if (auto it = matching.find(QCDOrder::LO);   it != matching.end()) return it->second;
    return {};
}

/**
 * @brief Build a @ref CustomCoefficientGroup from a lambda config.
 *
 * @param cfg User-facing group config.
 * @param adapters Adapter bundle used by the Wilson domain layer.
 * @return A ready-to-register group whose matching and running rules are wired.
 */
inline std::shared_ptr<CustomCoefficientGroup> make_custom_wilson_group(
    const CustomWilsonGroupConfig& cfg,
    WilsonGroupAdapterConfig adapters
) {
    const std::string name = cfg.display_name.empty() ? GroupMapper::str(cfg.group) : cfg.display_name;
    auto group = std::make_shared<CustomCoefficientGroup>(adapters, cfg.group, name, cfg.contribution);
    group->set_matching_storage_block(WilsonBlockNames::matching(cfg.group));
    group->set_wilson_type(cfg.contribution);

    for (const auto& coef_cfg : cfg.coefficients) {
        auto coef = std::make_shared<CustomWilson>(
            coef_cfg.id,
            WilsonBlockNames::matching(cfg.group),
            cfg.contribution
        );

        for (auto order : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
            auto it = coef_cfg.matching.find(order);
            if (it != coef_cfg.matching.end()) {
                coef->set_order_info(order, it->second.sources, it->second.compute, it->second.lhaid);
            } else {
                coef->set_order_info(
                    order,
                    {},
                    [](const ParamSrc&) { return 0.0; },
                    WCoefMapper::flha_full(coef_cfg.id, order, cfg.contribution)
                );
            }
        }

        group->add_coefficient(coef_cfg.id, coef);
    }

    if (!cfg.running.empty()) {
        for (const auto& [basis, per_order] : cfg.running) {
            group->add_sources(basis, per_order);
        }
    } else if (cfg.install_identity_running_if_empty) {
        std::map<QCDOrder, CoefficientGroupSources> per_order;
        for (auto order : {QCDOrder::LO, QCDOrder::NLO, QCDOrder::NNLO}) {
            CoefficientGroupSources s;
            s.sources = {{ParameterType::WILSON, {WilsonBlockNames::matching(cfg.group)}}};
            s.func = &custom_identity_running;
            per_order[order] = std::move(s);
        }
        group->add_sources(WilsonBasis::B_STANDARD, per_order);
    }

    return group;
}

#endif // CUSTOM_WILSON_LAMBDA_H
