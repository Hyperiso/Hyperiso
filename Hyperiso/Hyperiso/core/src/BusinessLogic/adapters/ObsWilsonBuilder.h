#ifndef OBS_WILSON_BUILDER_H
#define OBS_WILSON_BUILDER_H

#include "IObsWilsonBuilder.h"
#include "IWilsonBuilder.h"
#include "ObsWilsonProxy.h"
#include "WilsonProvider.h"
#include "Configs.h"
#include "Include.h"

/**
 * @file ObsWilsonBuilder.h
 * @brief Observable-side builder for Wilson coefficient services.
 *
 * This component is the "observable layer" entry-point to the Wilson pipeline.
 * It wraps a core @ref WilsonBuilder and exposes a simplified interface
 * compatible with the observable framework (taking an @ref AbstractConfig).
 *
 * Responsibilities:
 *  - Interpret an @ref AbstractConfig as a @ref WilsonBuildConfig.
 *  - Build Wilson coefficient groups if nothing exists yet, otherwise add them.
 *  - Provide an @ref IObsWilsonProxy used by observables to query coefficients.
 *
 * The build logic follows:
 *  - If the underlying core builder already owns a @ref CoefficientManager, we call
 *    @ref WilsonBuilder::add to append new groups to the existing manager.
 *  - Otherwise, we call @ref WilsonBuilder::build to create everything from scratch.
 *
 * @see WilsonBuilder
 * @see WilsonProvider
 * @see ObsWilsonProxy
 * @see IObsWilsonBuilder
 */
class ObsWilsonBuilder : public IObsWilsonBuilder {
public:
    /**
     * @brief Construct the observable builder using a core Wilson builder.
     * @param wil_builder Core Wilson builder responsible for group creation and configuration.
     *
     * The builder instance is stored and reused to either build or extend the
     * Wilson coefficient manager.
     */
    ObsWilsonBuilder(std::shared_ptr<WilsonBuilder> wil_builder) : wil_builder(std::move(wil_builder)) {}

    /**
     * @brief Build or extend the Wilson coefficient infrastructure from an observable config.
     *
     * The provided @p config must be a @ref WilsonBuildConfig (stored behind
     * an @ref AbstractConfig shared_ptr). It is downcast and then forwarded to the core builder.
     *
     * Behavior:
     *  - If a coefficient manager already exists, the config is appended via @ref WilsonBuilder::add.
     *  - Otherwise, the Wilson pipeline is created via @ref WilsonBuilder::build.
     *
     * @param config Observable configuration (expected to contain a WilsonBuildConfig).
     *
     * @warning This function assumes @p config can be safely cast to @ref WilsonBuildConfig.
     *          If the cast fails, dereferencing the result is undefined behavior.
     */
    void build(std::shared_ptr<AbstractConfig> config) override;

    /**
     * @brief Append a runtime/custom Wilson group to the underlying core builder.
     *
     * If the core builder has not created a coefficient manager yet, an empty
     * Wilson pipeline is first built using the scales/order from @p config.
     * This makes the method safe to call before any builtin Wilson group has
     * been requested by a decay.
     *
     * @param config Runtime/custom Wilson group configuration.
     */
    void add_custom_group(const CustomWilsonGroupConfig& config) override;

    /**
     * @brief Create and return a proxy object used by observables to query Wilson coefficients.
     *
     * The returned proxy is an @ref ObsWilsonProxy that wraps the core @ref WilsonProvider
     * obtained from the underlying @ref WilsonBuilder.
     *
     * @return A new @ref IObsWilsonProxy instance.
     */
    std::shared_ptr<IObsWilsonProxy> get_proxy() override;

private:
    /// Core Wilson builder (shared with the Wilson subsystem).
    std::shared_ptr<WilsonBuilder> wil_builder;
};

#endif