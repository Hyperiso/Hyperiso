#ifndef WILSONINTERFACE_H
#define WILSONINTERFACE_H

#include <map>

#include "Include.h"
#include "WilsonManager.h"
#include "AbstractConfig.h"
#include "WilsonBuilder.h"

/**
 * @brief User-facing API to build and query Wilson coefficients at matching and hadronic scales.
 *
 * @details
 * This is the *main* entry point for end users. It provides:
 *
 * 1) Lifecycle:
 *    - @ref build() to initialize the Wilson pipeline (groups, helpers, scales, blocks),
 *    - @ref addWilsonGroup() to extend an already built instance with additional groups,
 *    - @ref set_matching_scale() / @ref set_hadronic_scale() to update active scales.
 *
 * 2) Queries at MATCHING scale:
 *    - @ref getMatchingCoefficient() / @ref getM() :
 *        returns the coefficient at a given QCD order (LO/NLO/NNLO) for a contribution component.
 *    - @ref getFullMatchingCoefficient() / @ref getFM() :
 *        returns the coefficient including the perturbative sum up to the requested order.
 *
 * 3) Queries at HADRONIC (running) scale:
 *    - @ref getRunCoefficient() / @ref getR() :
 *        returns the evolved coefficient at a given QCD order (basis-dependent).
 *    - @ref getFullRunCoefficient() / @ref getFR() :
 *        returns the evolved coefficient including perturbative sum up to requested order.
 *
 * 4) Bulk convenience helpers:
 *    - per-order maps: @ref getSepOrderMatchingCoefficient() (@ref getSM),
 *                      @ref getSepOrderRunCoefficient() (@ref getSR)
 *    - full group maps: @ref getAllMatchingCoefficients() (@ref getAM),
 *                       @ref getAllRunCoefficients() (@ref getAR),
 *                       @ref getAllFullMatchingCoefficients() (@ref getAFM),
 *                       @ref getAllFullRunCoefficients() (@ref getAFR)
 *
 * Contribution semantics:
 * - @ref ContributionType::SM    : SM-only component.
 * - @ref ContributionType::BSM   : BSM-only component.
 * - @ref ContributionType::TOTAL : SM + BSM (composed by the manager depending on backend/input).
 *
 * Backend / MARTY caveat:
 * - When MARTY is enabled, the current pipeline computes Wilson coefficients at LO only.
 *   Requests at NLO/NNLO are expected to be downgraded to LO by the underlying provider.
 *   (The helper method ensure_mty_compat() exists, but the actual downgrade is enforced in WilsonProvider.)
 *
 * Threading / mutability:
 * - This interface mutates global-ish scale parameters through @ref ScaleSetter.
 *   If using multiple instances concurrently, ensure scale setters and the underlying parameter storage
 *   are safe for your execution model.
 *
 * Typical usage:
 * @code
 *   WilsonInterface W;
 *   WilsonBuildConfig cfg;
 *   cfg.groups = { GroupMapper::to_id(WGroup::B) };
 *   cfg.matching_scale = 160.0;
 *   cfg.hadronic_scale = 4.8;
 *   cfg.order = QCDOrder::NNLO;
 *   W.build(cfg);
 *
 *   auto C7_sm = W.getM(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::SM);
 *   auto C7_tot_full = W.getFM(WGroup::B, WCoef::C7, QCDOrder::NNLO, ContributionType::TOTAL);
 * @endcode
 *
 * @see WilsonBuilder
 * @see WilsonProvider
 * @see CoefficientManager
 * @see WilsonRequest
 */
class WilsonInterface {
private:
    std::shared_ptr<WilsonBuilder> builder;
    std::shared_ptr<WilsonProvider> provider;

    /// Patches registered before build(); applied immediately after build().
    std::vector<WilsonMatchingPatch> pending_matching_patches;

    /**
     * @brief Ensures requested QCD order is compatible with MARTY limitations.
     *
     * Current implementation of this helper logs and returns LO if MARTY is enabled
     * and a higher order is requested.
     *
     * Note: In the provided implementation, order downgrading is actually enforced
     * in @ref WilsonProvider::get() as well.
     */
    QCDOrder ensure_mty_compat(QCDOrder order);

    /// Whether build() has been called successfully.
    bool built {false};

public:
    WilsonInterface() = default;

    /**
     * @brief Builds the full Wilson pipeline and makes queries available.
     *
     * This constructs a @ref WilsonBuilder with the provided config and retrieves
     * the corresponding @ref WilsonProvider.
     *
     * After calling build(), all query methods become valid.
     */
    void build(WilsonBuildConfig config);

    /**
     * @brief Adds one or more Wilson groups to an existing built interface.
     *
     * Requires @ref build() to have been called first.
     * Delegates to @ref WilsonBuilder::add().
     */
    void addWilsonGroup(WilsonBuildConfig config);


    /**
     * @brief Register and initialize a runtime Wilson group described by lambdas.
     *
     * @details
     * This is the high-level entry point for user-defined Wilson coefficients.
     * The caller registers/resolves dynamic ids with @ref GroupMapper and
     * @ref WCoefMapper, fills a @ref CustomWilsonGroupConfig with matching
     * lambdas (and optional running lambdas), then passes it here.
     *
     * If the interface has not been built yet, an empty Wilson pipeline is
     * initialized first using the scales and order stored in @p config.
     *
     * @param config Runtime/custom Wilson group configuration.
     * @return Reference to `*this` for chaining.
     *
     * @see CustomWilsonGroupConfig
     * @see CustomWilsonCoefficientConfig
     */
    WilsonInterface& addCustomWilsonGroup(const CustomWilsonGroupConfig& config);

    /**
     * @brief Snake-case alias for @ref addCustomWilsonGroup.
     */
    WilsonInterface& add_custom_group(const CustomWilsonGroupConfig& config);

    /**
     * @brief Add an additive Wilson matching patch.
     *
     * If called before build(), the patch is queued and applied as soon as the
     * Wilson manager exists.  Otherwise it is applied immediately.
     */
    WilsonInterface& addMatchingPatch(const WilsonMatchingPatch& patch);

    /** @brief Snake-case alias for addMatchingPatch. */
    WilsonInterface& add_matching_patch(const WilsonMatchingPatch& patch);

    /** @brief Add several additive Wilson matching patches. */
    WilsonInterface& addMatchingPatches(const std::vector<WilsonMatchingPatch>& patches);

    /**
     * @brief Sets the matching scale (mu_W) in the global parameter system.
     *
     * This updates the value of the scale parameter corresponding to @ref ScaleType::MATCHING.
     * It does not automatically rebuild/combine blocks; dependent blocks will update on demand.
     */
    void set_matching_scale(double mu_W);

    /**
     * @brief Sets the hadronic scale (mu_h) in the global parameter system.
     *
     * This updates the value of the scale parameter corresponding to @ref ScaleType::HADRONIC.
     * It does not automatically rebuild/combine blocks; dependent blocks will update on demand.
     */
    void set_hadronic_scale(double mu_h);

    /**
     * @name Matching-scale single coefficient accessors
     *
     * The WGroup/WCoef overloads are kept for compatibility. The WGroupId/WCoefId
     * overloads are the preferred dynamic API for runtime/config-driven code.
     * @{
     */

    /**
     * @brief Returns the matching-scale coefficient at the requested QCD order.
     *
     * Dynamic-id overload. This avoids converting through the legacy static enum and
     * allows callers to resolve names with GroupMapper::id_of() and WCoefMapper::id_of().
     */
    scalar_t getMatchingCoefficient(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType cont_type);

    /// Backward-compatible enum overload for builtin groups and coefficients.
    scalar_t getMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type);

    /// Alias for @ref getMatchingCoefficient(WGroupId,WCoefId,QCDOrder,ContributionType).
    scalar_t getM(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType cont_type);

    /// Backward-compatible enum alias.
    scalar_t getM(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type);

    /**
     * @brief Returns the matching coefficient with perturbative summation up to `order`.
     *
     * Dynamic-id overload. The summation convention is implemented in the manager:
     *  C_full = Sum_{n=LO..order} C^{(n)} * (alpha_s(mu_W)/(4*pi))^(n-1).
     */
    scalar_t getFullMatchingCoefficient(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType cont_type);

    /// Backward-compatible enum overload for builtin groups and coefficients.
    scalar_t getFullMatchingCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type);

    /// Alias for @ref getFullMatchingCoefficient(WGroupId,WCoefId,QCDOrder,ContributionType).
    scalar_t getFM(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType cont_type);

    /// Backward-compatible enum alias.
    scalar_t getFM(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type);

    /** @} */

    /**
     * @name Hadronic-scale (running) single coefficient accessors
     * @{
     */

    /**
     * @brief Returns the hadronic/run coefficient at the requested order.
     *
     * Dynamic-id overload. @p basis selects the hadronic Wilson basis.
     */
    scalar_t getRunCoefficient(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Backward-compatible enum overload for builtin groups and coefficients.
    scalar_t getRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Alias for @ref getRunCoefficient(WGroupId,WCoefId,QCDOrder,ContributionType,WilsonBasis).
    scalar_t getR(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Backward-compatible enum alias.
    scalar_t getR(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /**
     * @brief Returns the full hadronic/run coefficient summed up to `order`.
     *
     * Dynamic-id overload. The summation convention is implemented in the manager:
     *  C_full = Sum_{n=LO..order} C^{(n)} * (alpha_s(mu_b)/(4*pi))^(n-1).
     */
    scalar_t getFullRunCoefficient(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Backward-compatible enum overload for builtin groups and coefficients.
    scalar_t getFullRunCoefficient(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Alias for @ref getFullRunCoefficient(WGroupId,WCoefId,QCDOrder,ContributionType,WilsonBasis).
    scalar_t getFR(WGroupId group, WCoefId coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Backward-compatible enum alias.
    scalar_t getFR(WGroup group, WCoef coeff, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /** @} */

    /**
     * @name Convenience: per-order maps for a fixed coefficient
     * @{
     */

    /**
     * @brief Returns a map {LO,NLO,NNLO} -> matching coefficient for a dynamic coefficient.
     */
    std::map<QCDOrder, scalar_t> getSepOrderMatchingCoefficient(WGroupId group, WCoefId coeff, ContributionType cont_type);

    /// Backward-compatible enum overload.
    std::map<QCDOrder, scalar_t> getSepOrderMatchingCoefficient(WGroup group, WCoef coeff, ContributionType cont_type);

    /// Alias for @ref getSepOrderMatchingCoefficient(WGroupId,WCoefId,ContributionType).
    std::map<QCDOrder, scalar_t> getSM(WGroupId group, WCoefId coeff, ContributionType cont_type);

    /// Backward-compatible enum alias.
    std::map<QCDOrder, scalar_t> getSM(WGroup group, WCoef coeff, ContributionType cont_type);

    /**
     * @brief Returns a map {LO,NLO,NNLO} -> hadronic coefficient for a dynamic coefficient.
     */
    std::map<QCDOrder, scalar_t> getSepOrderRunCoefficient(WGroupId group, WCoefId coeff, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Backward-compatible enum overload.
    std::map<QCDOrder, scalar_t> getSepOrderRunCoefficient(WGroup group, WCoef coeff, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Alias for @ref getSepOrderRunCoefficient(WGroupId,WCoefId,ContributionType,WilsonBasis).
    std::map<QCDOrder, scalar_t> getSR(WGroupId group, WCoefId coeff, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Backward-compatible enum alias.
    std::map<QCDOrder, scalar_t> getSR(WGroup group, WCoef coeff, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /** @} */

    /**
     * @name Convenience: full group maps (iterate all coefficients in the group)
     * @{
     */

    /// Returns {coef -> matching coefficient} for all coefficients in the group at a given order.
    std::map<WCoef, scalar_t> getAllMatchingCoefficients(WGroup group, QCDOrder order, ContributionType cont_type);

    /// Alias for @ref getAllMatchingCoefficients.
    std::map<WCoef, scalar_t> getAM(WGroup group, QCDOrder order, ContributionType cont_type);

    /// Returns {coef -> hadronic coefficient} for all coefficients in the group at a given order and basis.
    std::map<WCoef, scalar_t> getAllRunCoefficients(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Alias for @ref getAllRunCoefficients.
    std::map<WCoef, scalar_t> getAR(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Returns {coef -> full matching coefficient} (summed up to order) for all coefficients in the group.
    std::map<WCoef, scalar_t> getAllFullMatchingCoefficients(WGroup group, QCDOrder order, ContributionType cont_type);

    /// Alias for @ref getAllFullMatchingCoefficients.
    std::map<WCoef, scalar_t> getAFM(WGroup group, QCDOrder order, ContributionType cont_type);

    /// Returns {coef -> full hadronic coefficient} (summed up to order) for all coefficients in the group.
    std::map<WCoef, scalar_t> getAllFullRunCoefficients(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /// Alias for @ref getAllFullRunCoefficients.
    std::map<WCoef, scalar_t> getAFR(WGroup group, QCDOrder order, ContributionType cont_type, WilsonBasis basis=WilsonBasis::B_STANDARD);

    /** @} */
};

#endif // WILSONINTERFACE_H