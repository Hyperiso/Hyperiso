#ifndef OBSWILSONPROXY_H
#define OBSWILSONPROXY_H

#include "Include.h"
#include "IObsWilsonProxy.h"
#include "WilsonProvider.h"
#include "WilsonBuilder.h"

class ObsWilsonBuilder;

/**
 * @file ObsWilsonProxy.h
 * @brief Observable-facing proxy to query Wilson coefficients.
 *
 * This proxy is used by the observable layer to access Wilson coefficients without
 * exposing internal manager/build details. Internally it delegates to a
 * @ref WilsonProvider and uses @ref WilsonRequest as the transport structure.
 *
 * Supported queries:
 *  - Matching-scale coefficients (M): @ref getM
 *  - Matching-scale coefficients summed up to an order (FM): @ref getFM
 *  - Hadronic-scale (running) coefficients (R): @ref getR
 *  - Hadronic-scale coefficients summed up to an order (FR): @ref getFR
 *  - Convenience accessors returning maps across orders or all coefficients of a group:
 *    @ref getSM, @ref getSR, @ref getAM, @ref getAR, @ref getAFM, @ref getAFR.
 *
 * Basis handling:
 *  - A default basis (@ref WilsonBasis::B_STANDARD) is stored in the proxy.
 *  - The basis can be updated through @ref set_basis.
 *  - The basis is injected into every @ref WilsonRequest (even matching requests, for consistency).
 *
 * @note The proxy defaults the contribution to @ref ContributionType::TOTAL, which matches the
 *       typical observable need (“total” Wilson coefficient used in predictions).
 *
 * @see WilsonProvider
 * @see WilsonRequest
 * @see IObsWilsonProxy
 * @see ObsWilsonBuilder
 */
class ObsWilsonProxy : public IObsWilsonProxy {
public:
    /**
     * @brief Construct a proxy from an existing Wilson provider.
     * @param wilson_provider Provider used to resolve Wilson requests.
     */
    ObsWilsonProxy(std::shared_ptr<WilsonProvider> wilson_provider) : wil_p(wilson_provider) {}

    /**
     * @brief Get a matching-scale Wilson coefficient at a given QCD order (no order summation).
     *
     * This builds a @ref WilsonRequest with:
     *  - scale_type = @ref ScaleType::MATCHING
     *  - sum_qcd_orders = false
     *
     * @param group Wilson group enum (e.g. B, K, MesonMixing, ...).
     * @param coeff Wilson coefficient enum within the group.
     * @param order QCD order of the coefficient (LO/NLO/NNLO).
     * @param contribution Contribution type (SM/BSM/TOTAL). Default is TOTAL.
     * @return The complex Wilson coefficient value at matching scale.
     */
    complex_t getM(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) override;

    /**
     * @brief Get a matching-scale coefficient summed up to @p order (LO + ... + order).
     *
     * This sets sum_qcd_orders = true in the request. The underlying provider/manager
     * applies the appropriate perturbative expansion factorization.
     *
     * @param group Wilson group enum.
     * @param coeff Wilson coefficient enum.
     * @param order Maximum QCD order included in the sum.
     * @param contribution Contribution type (default TOTAL).
     * @return The complex coefficient summed up to the given order at matching scale.
     */
    complex_t getFM(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) override;

    /**
     * @brief Get a hadronic-scale (running) Wilson coefficient at a given QCD order (no summation).
     *
     * This builds a @ref WilsonRequest with:
     *  - scale_type = @ref ScaleType::HADRONIC
     *  - sum_qcd_orders = false
     *
     * The stored basis is used (see @ref set_basis).
     *
     * @param group Wilson group enum.
     * @param coeff Wilson coefficient enum.
     * @param order QCD order (LO/NLO/NNLO).
     * @param contribution Contribution type (default TOTAL).
     * @return The complex coefficient at hadronic scale.
     */
    complex_t getR(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) override;

    /**
     * @brief Get a hadronic-scale coefficient summed up to @p order (LO + ... + order).
     *
     * This sets sum_qcd_orders = true in the request and uses the stored basis.
     *
     * @param group Wilson group enum.
     * @param coeff Wilson coefficient enum.
     * @param order Maximum QCD order included in the sum.
     * @param contribution Contribution type (default TOTAL).
     * @return The complex coefficient summed up to the given order at hadronic scale.
     */
    complex_t getFR(WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) override;

    /**
     * @brief Return matching-scale coefficients for LO/NLO/NNLO separately.
     *
     * Convenience helper returning a map:
     *  - LO  -> getM(..., LO)
     *  - NLO -> getM(..., NLO)
     *  - NNLO-> getM(..., NNLO)
     *
     * @param group Wilson group enum.
     * @param coeff Wilson coefficient enum.
     * @param contribution Contribution type (default TOTAL).
     * @return Map keyed by QCD order.
     */
    std::map<QCDOrder, complex_t> getSM(WGroup group, WCoef coeff, ContributionType contribution=ContributionType::TOTAL);

    /**
     * @brief Return hadronic-scale coefficients for LO/NLO/NNLO separately.
     *
     * Convenience helper returning a map:
     *  - LO  -> getR(..., LO)
     *  - NLO -> getR(..., NLO)
     *  - NNLO-> getR(..., NNLO)
     *
     * @param group Wilson group enum.
     * @param coeff Wilson coefficient enum.
     * @param contribution Contribution type (default TOTAL).
     * @return Map keyed by QCD order.
     */
    std::map<QCDOrder, complex_t> getSR(WGroup group, WCoef coeff, ContributionType contribution=ContributionType::TOTAL);

    /**
     * @brief Return all matching-scale coefficients of a group at a given order.
     *
     * Iterates over @ref WCoefMapper::get_group(group) and queries @ref getM for each coefficient.
     *
     * @param group Wilson group enum.
     * @param order QCD order used for all coefficients.
     * @param contribution Contribution type (default TOTAL).
     * @return Map keyed by coefficient enum.
     */
    std::map<WCoef, complex_t> getAM(WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL);

    /**
     * @brief Return all hadronic-scale coefficients of a group at a given order.
     *
     * Iterates over @ref WCoefMapper::get_group(group) and queries @ref getR for each coefficient.
     *
     * @param group Wilson group enum.
     * @param order QCD order used for all coefficients.
     * @param contribution Contribution type (default TOTAL).
     * @return Map keyed by coefficient enum.
     */
    std::map<WCoef, complex_t> getAR(WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL);

    /**
     * @brief Return all full matching-scale coefficients of a group summed up to @p order.
     *
     * Iterates over @ref WCoefMapper::get_group(group) and queries @ref getFM for each coefficient.
     *
     * @param group Wilson group enum.
     * @param order Maximum QCD order in the sum for all coefficients.
     * @param contribution Contribution type (default TOTAL).
     * @return Map keyed by coefficient enum.
     */
    std::map<WCoef, complex_t> getAFM(WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL);

    /**
     * @brief Return all full hadronic-scale coefficients of a group summed up to @p order.
     *
     * Iterates over @ref WCoefMapper::get_group(group) and queries @ref getFR for each coefficient.
     *
     * @param group Wilson group enum.
     * @param order Maximum QCD order in the sum for all coefficients.
     * @param contribution Contribution type (default TOTAL).
     * @return Map keyed by coefficient enum.
     */
    std::map<WCoef, complex_t> getAFR(WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL);

    /**
     * @brief Recreate a builder interface usable from the observable layer.
     *
     * This retrieves the underlying core builder from the provider, casts it to
     * @ref WilsonBuilder, and wraps it into a new @ref ObsWilsonBuilder.
     *
     * @return A new observable builder sharing the same core builder.
     */
    std::shared_ptr<IObsWilsonBuilder> get_builder() override;

    /**
     * @brief Return the set of available bases for a given Wilson group id.
     *
     * Delegates to @ref WilsonProvider::get_bases.
     *
     * @param group Group identifier.
     * @return Set of supported bases.
     */
    std::unordered_set<WilsonBasis> get_bases(WGroupId) override;

    /**
     * @brief Set the default basis used for hadronic-scale queries.
     *
     * The stored basis is applied to every request (matching and hadronic for consistency),
     * but it is only physically relevant for hadronic blocks.
     *
     * @param basis New basis to use (e.g. B_STANDARD, SUSY basis, etc.).
     */
    void set_basis(WilsonBasis basis) override;

private:
    /// Core Wilson provider used to answer requests.
    std::shared_ptr<WilsonProvider> wil_p;

    /// Default basis applied to requests (mainly used at hadronic scale).
    WilsonBasis basis {WilsonBasis::B_STANDARD}; 
};

#endif // __OBSWILSONPROXY_H__
