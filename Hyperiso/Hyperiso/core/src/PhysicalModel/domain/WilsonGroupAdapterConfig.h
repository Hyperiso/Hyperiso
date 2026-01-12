#ifndef WILSON_GROUP_ADAPTER_CONFIG_H
#define WILSON_GROUP_ADAPTER_CONFIG_H

#include <memory>

#include "InterpretedParam.h"
#include "IBlockComposer.h"
#include "ICoreAPI.h"
#include "IMartyWilsonProxy.h"
#include "IParameterProxy.h"
#include "config.hpp"

/**
 * @file WilsonGroupAdapterConfig.h
 * @brief Adapter bundle wiring a CoefficientGroup to framework services.
 *
 * This header defines @ref WilsonGroupAdapterConfig, a small configuration struct
 * collecting all adapters and ports needed by the Wilson-domain layer.
 *
 * It is used by @ref CoefficientGroup to:
 * - read/write Wilson coefficients through a WILSON-scoped @ref IParameterProxy,
 * - register dependencies (dependent parameters / blocks) through an @ref IBlockComposer,
 * - query core framework state through @ref ICoreAPI (e.g. whether MARTY is used),
 * - optionally dispatch Wilson computations to a MARTY backend via @ref IMartyWilsonProxy.
 *
 * This indirection avoids hard-coding concrete implementations (HyperisoMaster, ParameterProxy,
 * WilsonParamComposer, MartyWilsonProxy, …) inside the domain objects.
 */

/**
 * @struct WilsonGroupAdapterConfig
 * @ingroup DomainModule
 * @brief Dependency injection bundle for CoefficientGroup.
 *
 * Required adapters:
 * - @ref wilson_proxy: provides access to coefficient values stored in the WILSON parameter space.
 * - @ref iblock_c: dependency composer used to register DependentParameters / blocks.
 * - @ref use_marty: boolean core API controlling whether MARTY is used.
 * - @ref marty_model_name: core API providing MARTY model name.
 * - @ref marty_model_path: core API providing MARTY model header path.
 *
 * Optional adapter:
 * - @ref marty_proxy: if provided, allows offloading Wilson coefficient computation/generation
 *   to MARTY (and tracking dependencies).
 *
 * Convenience:
 * - @ref sm_path points to the default SM model header under assets.
 */
struct WilsonGroupAdapterConfig {
    /**
     * @brief Constructs the adapter bundle.
     *
     * @param wilson_proxy     Proxy used to access WILSON values by (block, LhaID).
     * @param ibc             Block composer used to build dependent parameters/blocks.
     * @param use_marty       Core API returning whether MARTY is enabled (Model == MARTY).
     * @param marty_model_name Core API returning MARTY model name.
     * @param marty_model_path Core API returning MARTY model header path.
     * @param marty_proxy     Optional MARTY Wilson proxy (may be nullptr).
     */
    WilsonGroupAdapterConfig(std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy, std::shared_ptr<IBlockComposer> ibc,
    std::shared_ptr<ICoreAPI<bool>> use_marty, std::shared_ptr<ICoreAPI<std::string>> marty_model_name, std::shared_ptr<ICoreAPI<fs::path>> marty_model_path,
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr) 
    : wilson_proxy(wilson_proxy), iblock_c(ibc), use_marty(use_marty), marty_model_name(marty_model_name), marty_model_path(marty_model_path), marty_proxy(marty_proxy) {}

    /// Proxy used to read/write Wilson coefficients from Parameters (WILSON scope).
    std::shared_ptr<IParameterProxy<std::string, LhaID>> wilson_proxy;

     /// Composer used to register dependent parameters / dependent blocks.
    std::shared_ptr<IBlockComposer> iblock_c;

    /// Core API: true if MARTY is used for the current run.
    std::shared_ptr<ICoreAPI<bool>> use_marty;

    /// Core API: provides current MARTY model name.
    std::shared_ptr<ICoreAPI<std::string>> marty_model_name;

    /// Core API: provides current MARTY model header path.
    std::shared_ptr<ICoreAPI<fs::path>> marty_model_path;

    /// Optional MARTY Wilson proxy used to compute coefficients and fetch dependencies.
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy = nullptr;

    /**
     * @brief Default path to the SM MARTY header under assets.
     *
     * Used as a fallback/reference when building resolver/configuration.
     */
    fs::path sm_path = fs::path(std::string(project_assets_root.data())+"input_files/marty_model/sm.h");
};

#endif // WILSON_GROUP_ADAPTER_CONFIG_H