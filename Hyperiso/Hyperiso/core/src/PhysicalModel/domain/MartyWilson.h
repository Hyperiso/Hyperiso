#ifndef MARTY_WILSON_H
#define MARTY_WILSON_H

#include <complex>
#include <iostream>
#include <math.h>

#include "Wilson.h"
#include "config.hpp"
#include "DataFrame.h"
#include "CSVReader.h"
#include "InterpretedParam.h" //Only for template argument !!
#include "IMartyWilsonProxy.h"
#include "IMartyWilsonPathProxy.h"
#include "config.hpp"
#include "Utils.h"

/**
 * @file MartyWilson.h
 * @brief WilsonCoefficient implementation whose matching is produced by a MARTY backend.
 *
 * This header defines:
 * - @ref MartyWilsonConfig : configuration bundle describing how to obtain a coefficient
 *   from MARTY (model name/path, storage block, coefficient id, CSV location, proxy).
 * - @ref MartyWilson : a concrete @ref WilsonCoefficient whose LO matching value is
 *   computed by invoking a MARTY-based generator and reading the result from a CSV.
 *
 * High-level behavior
 * -------------------
 * MartyWilson is meant to plug into the same dependency/composition pipeline as
 * builtin coefficients:
 * - The framework composes dependent parameters for matching values (LO/NLO/NNLO LhaID).
 * - For MartyWilson, the LO compute functor:
 *   1) reads the current matching scale (EW_SCALE) from sources,
 *   2) asks @ref IMartyWilsonProxy to generate the coefficient for this scale,
 *   3) parses the CSV output (one row per Q_match),
 *   4) extracts the complex coefficient value for the requested coefficient name,
 *   5) lazily discovers and caches the list of parameter dependencies required by MARTY.
 *
 * Notes
 * -----
 * - Only LO provides a compute() function in this implementation; NLO/NNLO ids are
 *   still populated so the coefficient can be addressed consistently in blocks.
 * - Dependency discovery is done once (guarded by a size check on the ParamSrc raw map):
 *   the first time the functor is executed in the real dependency system, the proxy
 *   provides @ref IMartyWilsonProxy::get_dependencies() and @ref get_special_blocks().
 *
 * @see WilsonCoefficient
 * @see MatchingInfo
 * @see IMartyWilsonProxy
 * @see CoefficientRegistry (for factory registration)
 */

/**
 * @struct MartyWilsonConfig
 * @brief Configuration bundle for constructing a @ref MartyWilson coefficient.
 *
 * It describes:
 * - the MARTY model (name and header path),
 * - the output CSV location used by the MARTY generator,
 * - the coefficient LHA id (including order / contribution type encoded in parts),
 * - the destination block where matching values are stored,
 * - a proxy used to invoke the MARTY pipeline and query its dependencies.
 *
 * Default values target the SM model shipped in the assets directory.
 */
struct MartyWilsonConfig {
    /// Output model label used for generated libraries/CSV files (e.g. "SM", "THDM").
    std::string model_name{"SM"};

    /// C++ MARTY model class to instantiate. Usually equal to model_name, except
    /// for SM-like generation inside a BSM model.
    std::string generation_model_name{"SM"};

    /// If true, generated code instantiates generation_model_name but disables
    /// all particles not present in SM_Model.
    bool sm_like_filter{false};

    /// If true, generated MARTY code writes the BSM piece as TOTAL - SM-like.
    bool bsm_split_generation{false};
    
    /// Full coefficient id used to store the matching value (including order/type parts).
    LhaID coeff_id;

    /// Name of the block where matching coefficients are stored (e.g. "WC_B_MATCH").
    std::string storage_block;

    /// Path to the MARTY model header file.
    fs::path model_path{project_assets_root.data() + std::string() + "/input_files/marty_model/sm.h"};

    /// Proxy used to run MARTY and retrieve dependencies/special blocks.
    std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> marty_proxy;

    /// Path provider used to resolve MARTY generated outputs.
    std::shared_ptr<IMartyWilsonPathProxy> path_provider;

    /// Absolute path to the CSV output produced by MARTY for this model.
    std::string csv_path;

    /**
     * @brief Fallback path used only when legacy code constructs MartyWilsonConfig without a path provider.
     *
     * New code should pass @ref path_provider so the CSV path comes from
     * MemoryManager/APIPath::MARTY_TEMP_DIR rather than project_assets_root.
     */
    static fs::path fallback_csv_path(const std::string& model_name) {
        return fs::temp_directory_path() / "pyhyperiso" / "MartyTemp" / (model_name + "_wilson.csv");
    }

    /**
     * @brief Resolves and stores the CSV path using the path provider when available.
     */
    void refresh_csv_path() {
        if (path_provider) {
            csv_path = path_provider->wilson_csv_path(model_name).string();
        } else {
            csv_path = fallback_csv_path(model_name).string();
        }
    }

    /**
     * @brief Constructs config for a coefficient with a provided model path and proxy.
     *
     * @param id                 LHA identifier for the coefficient.
     * @param storage_block_name Destination block name where the value is stored.
     * @param model_path         Path to the MARTY model header file.
     * @param proxy              Proxy to run MARTY and query dependencies.
     * @param paths              Path provider resolving MARTY temp/CSV outputs.
     */
    MartyWilsonConfig(const LhaID& id,
                      const std::string& storage_block_name,
                      fs::path model_path,
                      std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> proxy,
                      std::shared_ptr<IMartyWilsonPathProxy> paths = nullptr)
        : coeff_id(id),
          storage_block(storage_block_name),
          model_path(model_path),
          marty_proxy(proxy),
          path_provider(paths) {
        refresh_csv_path();
    }

    /**
     * @brief Constructs config including an explicit model name.
     *
     * The CSV path is resolved from @ref IMartyWilsonPathProxy, typically to
     * `<cache>/MartyTemp/<model_name>_wilson.csv`.
     *
     * @param model_name         MARTY model class name.
     * @param id                 LHA identifier for the coefficient.
     * @param storage_block_name Destination block name where the value is stored.
     * @param model_path         Path to the MARTY model header file.
     * @param proxy              Proxy to run MARTY and query dependencies.
     * @param paths              Path provider resolving MARTY temp/CSV outputs.
     */
    MartyWilsonConfig(const std::string& model_name,
                      const LhaID& id,
                      const std::string& storage_block_name,
                      fs::path model_path,
                      std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> proxy,
                      std::shared_ptr<IMartyWilsonPathProxy> paths = nullptr)
        : model_name(model_name),
          generation_model_name(model_name),
          coeff_id(id),
          storage_block(storage_block_name),
          model_path(model_path),
          marty_proxy(proxy),
          path_provider(paths) {
        refresh_csv_path();
    }

    /**
     * @brief Constructs config with a separate output label and generated model.
     *
     * Used for SM-like contributions inside an extended MARTY model.
     */
    MartyWilsonConfig(const std::string& output_model_name,
                      const std::string& generation_model_name,
                      bool sm_like_filter,
                      bool bsm_split_generation,
                      const LhaID& id,
                      const std::string& storage_block_name,
                      fs::path model_path,
                      std::shared_ptr<IMartyWilsonProxy<InterpretedParam>> proxy,
                      std::shared_ptr<IMartyWilsonPathProxy> paths = nullptr)
        : model_name(output_model_name),
          generation_model_name(generation_model_name),
          sm_like_filter(sm_like_filter),
          bsm_split_generation(bsm_split_generation),
          coeff_id(id),
          storage_block(storage_block_name),
          model_path(model_path),
          marty_proxy(proxy),
          path_provider(paths) {
        refresh_csv_path();
    }
};

/**
 * @class MartyWilson
 * @ingroup DomainModule
 * @brief WilsonCoefficient whose matching value is computed through MARTY.
 *
 * MartyWilson derives from @ref WilsonCoefficient and installs a custom LO compute functor
 * into @ref WilsonCoefficient::matching_info:
 * - The compute functor calls @ref IMartyWilsonProxy::calculate() to generate values.
 * - It then reads the coefficient value from the configured CSV and returns it as scalar_t.
 * - It populates dependency sources on first meaningful execution by asking the proxy for:
 *   - @ref IMartyWilsonProxy::get_dependencies(),
 *   - @ref IMartyWilsonProxy::get_special_blocks().
 *
 * The contribution type is inferred from the LhaID "parts" (typically the 4th part).
 */
class MartyWilson : public WilsonCoefficient {
public:
    /**
     * @brief Constructs a Marty-backed Wilson coefficient.
     *
     * The base @ref WilsonCoefficient is initialized from @ref MartyWilsonConfig:
     * - coefficient name is inferred from coeff_id (via WCoefMapper),
     * - storage block is set from config.storage_block,
     * - contribution type is derived from coeff_id parts.
     *
     * The LO matching functor is installed and LO/NLO/NNLO LhaIDs are populated.
     *
     * @param config Marty configuration bundle.
     */
    MartyWilson(MartyWilsonConfig config);

    /**
     * @brief Returns the MARTY model name used by this coefficient.
     */
    std::string get_model() {
        return this->model;
    }

    /**
     * @brief Sets the MARTY model name used by this coefficient.
     */
    void set_model(std::string model) {
        this->model = model;
    }

    /**
     * @brief Polymorphic clone.
     *
     * @return A heap-allocated deep copy of this coefficient.
     */
    std::shared_ptr<WilsonCoefficient> clone() const override {
        return std::make_shared<MartyWilson>(*this);
    }

private:
    /// Name of the MARTY model to be used when generating the coefficient.
    std::string model{"SM"};

};

#endif // MARTY_WILSON_H