#ifndef CONFIGS_H
#define CONFIGS_H

#include "Include.h"

#include <utility>

/**
 * @file Configs.h
 * @brief Configuration structures used to control Wilson coefficient computations,
 *        strong coupling evaluations and mass calculations.
 *
 * This header defines a set of lightweight POD-like configuration types derived from
 * AbstractConfig. They encapsulate the parameters required by the numerical routines
 * (Wilson coefficient builders, α_s calculators, mass calculators, etc.).
 */

/**
 * @struct WilsonBuildConfig
 * @brief Configuration for building sets of Wilson coefficients.
 *
 * This structure defines how Wilson coefficients should be constructed:
 *   - which operator groups are included,
 *   - the matching scale between UV and effective theories,
 *   - the hadronic (low-energy) scale,
 *   - and the perturbative QCD order to be used.
 *
 * The configuration can be built either from a set of group identifiers
 * (WGroupId) or directly from a set of group enums (WGroup), which are
 * internally mapped to their corresponding identifiers.
 */
struct WilsonBuildConfig : public AbstractConfig {
    /**
     * @brief Set of Wilson operator group identifiers to be included in the build.
     */
    std::unordered_set<WGroupId> groups;

    /**
     * @brief Matching scale (in GeV) at which the high-energy theory
     *        is matched to the effective theory.
     */
    double matching_scale;

    /**
     * @brief Hadronic scale (in GeV) at which the Wilson coefficients
     *        are evaluated in the effective theory.
     */
    double hadronic_scale;

    /**
     * @brief Perturbative QCD order used for the evolution and matching of Wilson coefficients.
     *        Defaults to leading order (QCDOrder::LO).
     */
    QCDOrder order {QCDOrder::LO};

    /**
     * @brief Default constructor.
     *
     * Initializes an empty group set and leaves the scales uninitialized.
     * The QCD order is set to leading order (QCDOrder::LO).
     */
    WilsonBuildConfig() = default;

    /**
     * @brief Constructs a WilsonBuildConfig from a set of group identifiers.
     *
     * @param groups_         Set of Wilson group identifiers to include.
     * @param matching_scale_ Matching scale at which matching is performed.
     * @param hadronic_scale_ Hadronic scale at which coefficients are evaluated.
     * @param order_          Perturbative QCD order (default: QCDOrder::LO).
     */
    WilsonBuildConfig(
        std::unordered_set<WGroupId> groups_,
        double matching_scale_,
        double hadronic_scale_,
        QCDOrder order_ = QCDOrder::LO
    ) :
        groups(std::move(groups_)),
        matching_scale(matching_scale_),
        hadronic_scale(hadronic_scale_),
        order(order_) {}

    /**
     * @brief Constructs a WilsonBuildConfig from a set of group enums.
     *
     * The provided WGroup values are mapped internally to their corresponding
     * WGroupId using GroupMapper::to_id().
     *
     * @param groups_         Set of Wilson groups to include.
     * @param matching_scale_ Matching scale at which matching is performed.
     * @param hadronic_scale_ Hadronic scale at which coefficients are evaluated.
     * @param order_          Perturbative QCD order (default: QCDOrder::LO).
     */
    WilsonBuildConfig(
        std::unordered_set<WGroup> groups_,
        double matching_scale_,
        double hadronic_scale_,
        QCDOrder order_ = QCDOrder::LO
    ) :
        matching_scale(matching_scale_),
        hadronic_scale(hadronic_scale_),
        order(order_) { for (auto group : groups_) { groups.insert(GroupMapper::to_id(group));}}
};

/**
 * @struct WilsonRequest
 * @brief Dynamic configuration for requesting a specific Wilson coefficient.
 *
 * Wilson requests are stored with the dynamic identifier layer:
 *   - @ref WGroupId for the Wilson group,
 *   - @ref WCoefId for the Wilson coefficient.
 *
 * This keeps the provider compatible with both:
 *   - builtin enum calls, through the compatibility constructor taking
 *     @ref WGroup and @ref WCoef,
 *   - runtime/config/CLI-driven calls, through the constructor taking
 *     @ref WGroupId and @ref WCoefId.
 *
 * The manager underneath already consumes string names, so using ids here avoids
 * converting through the static enum when the caller has a dynamic symbol.
 */
struct WilsonRequest : public AbstractConfig {
    /**
     * @brief Wilson operator group identifier associated with the request.
     *
     * For builtin groups, this is usually obtained with GroupMapper::to_id(WGroup::B).
     * For dynamic/config-driven code, prefer GroupMapper::id_of("B") or any registered alias.
     */
    WGroupId group;

    /**
     * @brief Wilson coefficient identifier associated with the request.
     *
     * For builtin coefficients, this is usually obtained with WCoefMapper::to_id(WCoef::C7).
     * For dynamic/config-driven code, prefer WCoefMapper::id_of("C7") or any registered alias.
     */
    WCoefId coefficient;

    /**
     * @brief Perturbative QCD order at which the coefficient is evaluated.
     */
    QCDOrder order {QCDOrder::LO};

    /**
     * @brief Type of contribution requested, e.g. SM, BSM or TOTAL.
     */
    ContributionType contribution {ContributionType::TOTAL};

    /**
     * @brief Scale at which the Wilson coefficient is evaluated.
     */
    ScaleType scale_type {ScaleType::HADRONIC};

    /**
     * @brief If true, contributions from different QCD orders are summed up to @ref order.
     */
    bool sum_qcd_orders {false};

    /**
     * @brief Optional Wilson basis used for hadronic-scale coefficients.
     */
    std::optional<WilsonBasis> basis;

    /**
     * @brief Constructs a WilsonRequest from dynamic identifiers.
     *
     * Use this overload for runtime/config driven code where the group and coefficient
     * have already been resolved by the mapper layer.
     *
     * @param group Wilson group id.
     * @param coefficient Wilson coefficient id.
     * @param order Perturbative QCD order.
     * @param contribution Contribution type to request.
     * @param scale_type Matching or hadronic scale.
     * @param sum_qcd_orders Whether to sum QCD orders up to @p order.
     */
    WilsonRequest(WGroupId group, WCoefId coefficient, QCDOrder order,
                  ContributionType contribution, ScaleType scale_type, bool sum_qcd_orders) :
        group(std::move(group)), coefficient(std::move(coefficient)), order(order),
        contribution(contribution), scale_type(scale_type), sum_qcd_orders(sum_qcd_orders) {}

    /**
     * @brief Backward-compatible constructor from builtin enums.
     *
     * This preserves existing calls such as:
     * @code
     * WilsonRequest{WGroup::B, WCoef::C7, QCDOrder::NLO,
     *                ContributionType::TOTAL, ScaleType::HADRONIC, true};
     * @endcode
     * while storing the request internally as dynamic ids.
     */
    WilsonRequest(WGroup group, WCoef coefficient, QCDOrder order,
                  ContributionType contribution, ScaleType scale_type, bool sum_qcd_orders) :
        WilsonRequest(GroupMapper::to_id(group), WCoefMapper::to_id(coefficient),
                      order, contribution, scale_type, sum_qcd_orders) {}
};

/**
 * @struct AlphasConfig
 * @brief Configuration for evaluating the strong coupling constant \f$\alpha_s\f$.
 *
 * This structure specifies:
 *   - the renormalization scale,
 *   - the mass scheme used for the bottom and top quarks.
 */
struct AlphasConfig : public AbstractConfig {
    /**
     * @brief Renormalization scale at which \f$\alpha_s\f$ is evaluated.
     */
    double scale;

    /**
     * @brief Mass scheme used for the bottom-quark mass (e.g. pole, MSbar).
     */
    MassType m_b_type;

    /**
     * @brief Mass scheme used for the top-quark mass (e.g. pole, MSbar).
     */
    MassType m_t_type;

    /**
     * @brief Constructs an AlphasConfig with the given parameters.
     *
     * @param scale      Renormalization scale for \f$\alpha_s\f$.
     * @param m_b_type   Bottom-quark mass scheme.
     * @param m_t_type   Top-quark mass scheme.
     */
    AlphasConfig(double scale, MassType m_b_type, MassType m_t_type) : scale(scale), m_b_type(m_b_type), m_t_type(m_t_type) {}
};

/**
 * @struct MassConfig
 * @brief Configuration for computing a particle mass at a given scale.
 *
 * This configuration extends AlphasConfig by adding:
 *   - the PDG identifier of the particle whose mass is evaluated.
 *
 * It is typically used to request running or scheme-converted masses
 * for specific particles at the given scale.
 */
struct MassConfig : public AlphasConfig {
    /**
     * @brief PDG identifier of the particle whose mass is requested.
     */
    int pdg_id;

    /**
     * @brief Constructs a MassConfig with the given parameters.
     *
     * @param pdg_id     PDG identifier of the particle.
     * @param scale      Renormalization scale at which the mass is evaluated.
     * @param m_b_type   Bottom-quark mass scheme (passed to AlphasConfig).
     * @param m_t_type   Top-quark mass scheme (passed to AlphasConfig).
     */
    MassConfig(int pdg_id, double scale, MassType m_b_type, MassType m_t_type) : AlphasConfig(scale, m_b_type, m_t_type), pdg_id(pdg_id) {}
};

#endif // CONFIGS_H
