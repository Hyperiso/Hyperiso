#ifndef CUSTOM_WILSON_H
#define CUSTOM_WILSON_H

#include <functional>
#include <memory>
#include <string>
#include <unordered_set>

#include "Wilson.h"

/**
 * @file CustomWilson.h
 * @brief User-definable WilsonCoefficient with explicit per-order matching definitions.
 *
 * This header defines @ref CustomWilson, a convenience concrete implementation of
 * @ref WilsonCoefficient intended for:
 * - quick prototyping,
 * - external/custom EFT extensions,
 * - unit tests,
 * where the user provides the dependency sources and compute functor directly.
 *
 * Unlike built-in or MARTY-backed coefficients, CustomWilson does not assume any
 * physics computation by itself. Instead, clients fill @ref MatchingInfo per QCD order
 * using @ref set_order_info().
 *
 * Typical usage
 * -------------
 * @code
 * auto c = std::make_shared<CustomWilson>(LhaID(1,2,0,0), "WCUSTOM");
 * c->set_order_info(QCDOrder::LO,
 *                   { ParamId{ParameterType::SM,"SMINPUTS",1} },
 *                   [](const ParamSrc& src){ return src.get_val("SMINPUTS",1); },
 *                   LhaID(1,2,0,0));
 * @endcode
 *
 * @see WilsonCoefficient
 * @see MatchingInfo
 */

/**
 * @class CustomWilson
 * @ingroup DomainModule
 * @brief Concrete WilsonCoefficient whose matching rules are provided by the user.
 *
 * The constructor initializes the base @ref WilsonCoefficient with an explicit LhaID,
 * storage block and contribution type, then resets per-order matching_info entries.
 * Users can then provide (sources, compute, lhaid) for each desired @ref QCDOrder via
 * @ref set_order_info().
 */
class CustomWilson : public WilsonCoefficient {
public:
    /**
     * @brief Constructs a custom Wilson coefficient.
     *
     * @param name          Base identifier for the coefficient (full LHA id).
     * @param storage_block Block where the coefficient is stored.
     * @param type          Contribution type (SM/BSM/TOTAL).
     */
    CustomWilson(const LhaID& name,
                 const std::string& storage_block,
                 ContributionType type = ContributionType::SM)
    : WilsonCoefficient(name, storage_block, type)
    {
        std::cout << "Creating CustomWilson: " << name << std::endl;
        matching_info[QCDOrder::LO] = MatchingInfo();
        matching_info[QCDOrder::NLO] = MatchingInfo();
        matching_info[QCDOrder::NNLO] = MatchingInfo();

        // this->max_order = max_order;
        this->type = type;
    }

    /**
     * @brief Constructs a custom coefficient from a dynamic Wilson coefficient id.
     *
     * This overload is the preferred entry point for runtime/user-defined Wilson
     * coefficients. The coefficient name is set to @ref WCoefMapper::str(id), so
     * it can be queried later through @ref WCoefId without requiring a legacy
     * @ref WCoef enum value. Per-order LHA ids are initialized with
     * @ref WCoefMapper::flha_full(WCoefId,QCDOrder,ContributionType).
     *
     * @param id            Dynamic coefficient id. It must have an FLHA base key
     *                      registered in @ref WCoefMapper.
     * @param storage_block Block where matching values will be composed.
     * @param type          Contribution slot computed by this coefficient.
     */
    CustomWilson(WCoefId id,
                 const std::string& storage_block,
                 ContributionType type = ContributionType::SM)
    : WilsonCoefficient(WCoefMapper::flha_full(id, QCDOrder::LO, type), storage_block, type)
    {
        this->coeffName = WCoefMapper::str(id);
        this->type = type;

        matching_info[QCDOrder::LO]   = MatchingInfo(WCoefMapper::flha_full(id, QCDOrder::LO, type));
        matching_info[QCDOrder::NLO]  = MatchingInfo(WCoefMapper::flha_full(id, QCDOrder::NLO, type));
        matching_info[QCDOrder::NNLO] = MatchingInfo(WCoefMapper::flha_full(id, QCDOrder::NNLO, type));
    }

    /**
     * @brief Defines sources and compute functor for a specific QCD order.
     *
     * @param order    Order to configure (LO/NLO/NNLO).
     * @param sources  Set of ParamId dependencies required for this computation.
     * @param compute  Function mapping a ParamSrc to the coefficient value.
     * @param lhaid    LhaID used to identify the stored parameter for this order.
     * @return Reference to *this (fluent API).
     */
    CustomWilson& set_order_info(
        QCDOrder order,
        std::unordered_set<ParamId> sources,
        std::function<scalar_t(const ParamSrc&)> compute,
        LhaID lhaid
    ) {
        matching_info[order] = MatchingInfo(std::move(sources), std::move(compute), std::move(lhaid));
        // if (order > max_order) max_order = order;
        return *this;
    }

    /**
     * @brief Sets the contribution type (fluent API).
     */
    CustomWilson& with_type(ContributionType t) { this->type = t; return *this; }

    /**
     * @brief Sets the storage block name (fluent API).
     */
    CustomWilson& with_storage_block(std::string block) { this->storage_block = std::move(block); return *this; }

    /**
     * @brief Polymorphic clone.
     *
     * @return A heap-allocated deep copy of this coefficient.
     */
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<CustomWilson>(*this);
    }
};

#endif // CUSTOM_WILSON_H
