#ifndef WILSONBUILDER_H
#define WILSONBUILDER_H

#include <memory>

#include "IWilsonBuilder.h"
#include "WilsonProvider.h"
#include "Include.h"
#include "WilsonManager.h"
#include "Configs.h"
#include "WilsonParamComposer.h"
#include "ScaleSetter.h"
#include "WilsonParametersHelper.h"
#include "ModelAPI.h"
#include "UseMarty.h"
#include "HasWilsonAPI.h"
#include "SUSYParametersHelper.h"
#include "THDMParametersHelper.h"
#include "MartyWilsonProxy.h"
#include "MartyModelNameAPI.h"
#include "MartyModelPathAPI.h"
#include "MartyWilsonPathProxy.h"
#include "SMFromHypProxy.h"
#include "CustomWilsonLambda.h"

/**
 * @brief High-level factory assembling a complete Wilson computation pipeline.
 *
 * WilsonBuilder is responsible for wiring together:
 * - The parameter composition engine (@ref WilsonParamComposer / @ref IBlockComposer),
 * - Wilson parameter helper blocks per model (SM always, plus THDM/SUSY when selected),
 * - The @ref CoefficientRegistry (factories for individual Wilson coefficients),
 * - The @ref CoefficientGroupBuilder (builds groups from @ref GroupDefinitions),
 * - The @ref CoefficientManager (registers groups, composes matching/hadronic blocks),
 * - Marty integration (optional proxy, model path/name APIs),
 * - Runtime APIs: model selection, FWCOEF presence, scale setters.
 *
 * Typical usage:
 * @code
 *   WilsonBuildConfig cfg;
 *   cfg.groups = { GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::K) };
 *   cfg.matching_scale = 160.0;
 *   cfg.hadronic_scale = 4.8;
 *   cfg.order = QCDOrder::NLO;
 *
 *   auto builder  = std::make_shared<WilsonBuilder>(cfg);
 *   auto provider = builder->get_wilson_provider();
 * @endcode
 *
 * Notes on Marty:
 * - When Marty is enabled, the builder currently forces coefficients to LO only
 *   (and will warn if a higher order is requested).
 * - Contribution type differs depending on backend: Marty tends to provide TOTAL,
 *   while builtin flow often separates SM vs BSM depending on selected model.
 *
 * @see CoefficientRegistry
 * @see CoefficientGroupBuilder
 * @see CoefficientManager
 * @see WilsonProvider
 */
class WilsonBuilder : public IWilsonBuilder<WilsonBuildConfig, WilsonProvider>, public std::enable_shared_from_this<WilsonBuilder> { 
public:
    /**
     * @brief Constructs and immediately builds the Wilson pipeline from a config.
     */
    WilsonBuilder(WilsonBuildConfig config);

    /**
     * @brief Constructs a builder from an existing coefficient manager (advanced usage).
     */
    WilsonBuilder(std::shared_ptr<CoefficientManager> manager);

    /**
     * @brief Build a new Wilson pipeline from scratch.
     *
     * Steps performed (as per implementation):
     * - create a shared block composer,
     * - initialize SM helper blocks for each requested group,
     * - depending on current model, initialize THDM/SUSY helper blocks (if applicable),
     * - create parameter proxies (Wilson + SM),
     * - set up runtime flags / APIs (use Marty? has Wilson input? model API),
     * - create Marty proxy if needed,
     * - create adapters and coefficient registry, build groups requested in config,
     * - enforce Marty LO-only if needed,
     * - build and initialize the CoefficientManager (matching + hadronic blocks).
     */
    void build(WilsonBuildConfig config) override;

    /**
     * @brief Add additional Wilson groups to an existing manager.
     *
     * This method:
     * - (re)initializes helper blocks for the new groups,
     * - builds the new group(s),
     * - registers them into the existing manager,
     * - initializes matching and hadronic blocks for each new group.
     *
     * Notes:
     * - This does not rebuild the whole pipeline; it extends the existing manager.
     * - Scales are updated before composing new blocks.
     */
    void add(WilsonBuildConfig config) override;

    /**
     * @brief Register and initialize a user-defined Wilson group built from lambdas.
     *
     * The config provides dynamic ids, matching lambdas for coefficients and optional
     * running lambdas for the group. The builder turns it into a
     * @ref CustomCoefficientGroup, registers it in the existing
     * @ref CoefficientManager, then initializes matching and hadronic blocks at the
     * scales/order stored in the config.
     *
     * @param config Runtime/custom Wilson group configuration.
     */
    void add_custom_group(const CustomWilsonGroupConfig& config);

    /**
     * @brief Add one additive matching patch to the active Wilson manager.
     */
    void add_matching_patch(const WilsonMatchingPatch& patch);

    /**
     * @brief Add several additive matching patches to the active Wilson manager.
     */
    void add_matching_patches(const std::vector<WilsonMatchingPatch>& patches);

    /**
     * @brief Returns a provider facade around this builder.
     *
     * The provider offers runtime coefficient access using a request/config object.
     */
    std::shared_ptr<WilsonProvider> get_wilson_provider();

    /**
     * @brief Returns the underlying coefficient manager.
     */
    std::shared_ptr<CoefficientManager> get_coefficient_manager();

private:
    /**
     * @brief Model-specific Wilson parameter helpers (SM always present).
     *
     * Helpers are responsible for composing auxiliary blocks (WPARAM_*, evolution matrices...)
     * via the shared block composer.
     */
    std::map<Model, std::shared_ptr<IWilsonParameterHelper>> wilson_param_helpers;

    /// Active coefficient manager (owns groups and composed blocks).
    std::shared_ptr<CoefficientManager> cm;

    /**
     * @brief Last adapter bundle used to build coefficient groups.
     *
     * It is cached so runtime custom groups can be assembled with the same proxies,
     * composer and MARTY hooks as builtin groups.
     */
    std::shared_ptr<WilsonGroupAdapterConfig> current_group_adapters;

    /// Cached MARTY path proxy matching current_group_adapters; reused by add().
    std::shared_ptr<IMartyWilsonPathProxy> current_marty_paths;
};

#endif // IWILSONBUILDER_H
