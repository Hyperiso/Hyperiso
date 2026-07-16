#ifndef WILSONPROVIDER_H
#define WILSONPROVIDER_H

#include <memory>

#include "IWilsonProvider.h"
#include "Include.h"
#include "WilsonManager.h"
#include "Configs.h"


class WilsonBuilder;

/**
 * @brief Thin query layer: converts user requests into CoefficientManager calls.
 *
 * @details
 * WilsonProvider is an execution/query façade:
 * - It keeps a reference to the @ref WilsonBuilder and the built @ref CoefficientManager.
 * - It exposes a single entry point `get(AbstractConfig)` used by the framework interfaces.
 *
 * The provider interprets @ref WilsonRequest:
 * - selects scale type (matching or hadronic),
 * - optionally sums QCD orders up to the requested order,
 * - selects basis for hadronic coefficients.
 *
 * MARTY compatibility:
 * - If MARTY is enabled and an order > LO is requested, the provider logs a warning
 *   and downgrades the request to LO (as per current backend limitation/policy).
 *
 * @see WilsonInterface
 * @see WilsonBuilder
 * @see CoefficientManager
 * @see WilsonRequest
 */
class WilsonProvider : public IWilsonProvider<WilsonBuilder> {
public:
    /**
     * @brief Constructs a provider from an existing builder.
     *
     * @param builder Non-null WilsonBuilder. The provider will also obtain the
     *                builder's coefficient manager.
     */
    WilsonProvider(std::shared_ptr<WilsonBuilder> builder);

    /**
     * @brief Returns the requested Wilson coefficient as scalar_t (complex).
     *
     * The @p config is expected to be a WilsonRequest wrapped as AbstractConfig.
     * Depending on request fields:
     * - scale_type == MATCHING  -> getMatchingCoefficient or getFullMatchingCoefficient
     * - scale_type == HADRONIC  -> getRunCoefficient or getFullRunCoefficient
     *
     * If Marty is enabled and request.order > LO, access is downgraded to LO with a warning.
     *
     * @param config Abstract configuration (must contain a WilsonRequest).
     * @return The requested coefficient value.
     */
    scalar_t get(std::shared_ptr<AbstractConfig> config) override;

    /**
     * @brief Returns the underlying builder.
     */
    std::shared_ptr<WilsonBuilder> get_builder() override;

    /**
     * @brief Returns the set of bases supported by a Wilson group.
     */
    std::unordered_set<WilsonBasis> get_bases(WGroupId group) override;

private:
    std::shared_ptr<WilsonBuilder> builder;
    std::shared_ptr<CoefficientManager> cm;
};

#endif // WILSONPROVIDER_H
