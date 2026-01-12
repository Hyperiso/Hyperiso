#ifndef IOBSWILSONPROXY_H
#define IOBSWILSONPROXY_H

#include "Include.h"

class IObsWilsonBuilder;

/**
 * @file IObsWilsonProxy.h
 * @brief Interface for querying Wilson coefficients from the observable layer.
 *
 * This interface is designed for observable computations. It exposes:
 *  - matching-scale accessors (M / FM),
 *  - hadronic-scale (running) accessors (R / FR),
 *  - bulk helpers returning maps over QCD orders or over all coefficients of a group,
 *  - access to the associated builder and supported bases.
 *
 * Naming conventions:
 *  - M  : matching-scale coefficient at a single QCD order.
 *  - FM : matching-scale coefficient summed up to a QCD order.
 *  - R  : hadronic-scale (running) coefficient at a single QCD order.
 *  - FR : hadronic-scale coefficient summed up to a QCD order.
 *  - SM/SR : "separated by order" maps for M/R.
 *  - AM/AR/AFM/AFR : "all coefficients" maps for M/R/FM/FR.
 *
 * @see ObsWilsonProxy
 * @see IObsWilsonBuilder
 */
class IObsWilsonProxy {
public:
    virtual ~IObsWilsonProxy() = default;

    /// @name Single coefficient accessors
    ///@{

    /**
     * @brief Matching-scale coefficient at a given order (no summation).
     */
    virtual complex_t                       getM    (WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) = 0;

    /**
     * @brief Matching-scale coefficient summed up to the given order.
     */
    virtual complex_t                       getFM   (WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) = 0;

    /**
     * @brief Hadronic-scale coefficient at a given order (no summation).
     */
    virtual complex_t                       getR    (WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) = 0;

    /**
     * @brief Hadronic-scale coefficient summed up to the given order.
     */
    virtual complex_t                       getFR   (WGroup group, WCoef coeff, QCDOrder order, ContributionType contribution=ContributionType::TOTAL) = 0;

    ///@}

    /// @name Convenience helpers (maps)
    ///@{

    /**
     * @brief Matching coefficients for LO/NLO/NNLO returned as a map keyed by QCD order.
     */
    virtual std::map<QCDOrder, complex_t>   getSM   (WGroup group, WCoef coeff, ContributionType contribution=ContributionType::TOTAL)                 = 0;

    /**
     * @brief Hadronic coefficients for LO/NLO/NNLO returned as a map keyed by QCD order.
     */
    virtual std::map<QCDOrder, complex_t>   getSR   (WGroup group, WCoef coeff, ContributionType contribution=ContributionType::TOTAL)                 = 0;

    /**
     * @brief All matching coefficients for a given group and order, keyed by coefficient enum.
     */
    virtual std::map<WCoef, complex_t>      getAM   (WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL)              = 0;

    /**
     * @brief All hadronic coefficients for a given group and order, keyed by coefficient enum.
     */
    virtual std::map<WCoef, complex_t>      getAR   (WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL)              = 0;

    /**
     * @brief All full matching coefficients for a given group (summed up to @p order).
     */
    virtual std::map<WCoef, complex_t>      getAFM  (WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL)              = 0;

    /**
     * @brief All full hadronic coefficients for a given group (summed up to @p order).
     */
    virtual std::map<WCoef, complex_t>      getAFR  (WGroup group, QCDOrder order, ContributionType contribution=ContributionType::TOTAL)              = 0;

    ///@}

    /// @name Builder / basis management
    ///@{

    /**
     * @brief Get a builder instance compatible with the observable layer.
     *
     * Intended for advanced workflows where observables want to trigger a (re)build or add groups.
     */
    virtual std::shared_ptr<IObsWilsonBuilder> get_builder() = 0;

    /**
     * @brief List available bases for a given Wilson group id.
     */
    virtual std::unordered_set<WilsonBasis> get_bases(WGroupId) = 0;

    /**
     * @brief Set the default basis used for hadronic-scale requests.
     */
    virtual void set_basis(WilsonBasis basis) = 0;

    ///@}
};

#endif // IOBSWILSONPROXY_H
