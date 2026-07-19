#ifndef MARTY_INTERFACE_H
#define MARTY_INTERFACE_H

#include <memory>
#include <filesystem>
#include <iostream>
#include <algorithm>
#include <set>
#include <unordered_set>
#include <unordered_map>
#include <optional>

#include "config.hpp"
#include "FileNameManager.h"
#include "GeneralModelModifier.h"
#include "CodeGenerator.h"
#include "GppCompilerStrategy.h"
#include "MakeCompilerStrategy.h"
#include "Logger.h"
#include "GeneralNumModelModifier.h"
#include "ICoreAPI.h"
#include "IMartyParameterProxy.h"

/**
 * @file MartyInterface.h
 * @brief High-level façade around MARTY-based code generation and compilation.
 *
 * This header defines ::MartyInterface, a convenience class that:
 *  - generates model-dependent C++ code from MARTY templates,
 *  - compiles and runs the generated code,
 *  - builds numeric libraries for Wilson coefficient evaluation,
 *  - exposes the set of parameter dependencies discovered during generation.
 *
 * It bundles together template management, mapping resolution, and
 * compiler strategies into a single entry point for the rest of the framework.
 */

/**
 * @class MartyInterface
 * @ingroup CodeGenerationModule
 * @brief High-level interface to MARTY-based code generation and execution.
 *
 * ::MartyInterface orchestrates the following steps:
 *  - non-numeric template generation using ::GeneralModelModifier and
 *    ::NonNumericTemplateManager,
 *  - compilation and execution of the generated code through
 *    ::GppCompilerStrategy,
 *  - numeric library generation (parameters, CSV helper, etc.) through
 *    ::GeneralNumModelModifier and ::NumericTemplateManager,
 *  - compilation and execution of numeric libraries with
 *    ::MakeCompilerStrategy,
 *  - bookkeeping of parameter dependencies (see ::get_dependencies()).
 *
 * It can be constructed either with explicit dependencies (core API,
 * parameter proxies, interpreter ports) or using the default constructor,
 * which wires default implementations together.
 */
class MartyInterface {
public:
    /**
     * @brief Constructs a MartyInterface with explicit dependencies.
     *
     * @param core_api  Shared pointer to the core API exposing the active ::Model.
     * @param sm_proxy  Proxy used to access SM parameters (blocks, LHA IDs).
     * @param ports     Factory used to build interpreter ports (resolvers, mapping adapters).
     * @param bsm_proxy Optional proxy for BSM parameters; may be null if not needed.
     */
    MartyInterface(std::shared_ptr<ICoreAPI<Model>> core_api, std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> sm_proxy, std::shared_ptr<IInterpreterPortsFactory> ports, std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> bsm_proxy = nullptr) : core_api(core_api), param_proxy_sm(sm_proxy), ports(ports), param_proxy_bsm(bsm_proxy)
     {}
    
    /**
     * @brief Default constructor wiring default implementations.
     *
     * Internally:
     *  - builds a ::ModelAPI for model access,
     *  - creates ::MartyParameterProxy instances for SM and BSM,
     *  - uses ::DefaultInterpreterPortsFactory for mapping resolution.
     */
    MartyInterface();

    /**
     * @brief Generates a non-numeric MARTY-powered C++ template.
     *
     * This creates a C++ file (e.g. for calling a Wilson calculator) by:
     *  - configuring a ::GeneralModelModifier,
     *  - using a ::NonNumericTemplateManager and ::CodeGenerator.
     *
     * @param wilson     Name of the Wilson basis.
     * @param model      Target model name (e.g. "SM", "THDM").
     * @param model_path Filesystem path to the model header used in the template.
     */
    void generate(std::string wilson, std::string model, std::string model_path);

    /**
     * @brief Generates code with a separate output label and target model.
     *
     * This is used for SM-like contributions in an extended model: output files
     * and libraries keep the ``SM`` label, but the generated code instantiates
     * the target BSM model and may filter out non-SM particles.
     */
    void generate(std::string wilson,
                  std::string output_model,
                  std::string target_model,
                  std::string model_path,
                  bool sm_like_filter,
                  bool bsm_split_generation = false,
                  bool full_target_generation = false);

    /**
     * @brief Compiles and runs the non-numeric generated code.
     *
     * Uses ::GppCompilerStrategy to:
     *  - compile the generated C++ file,
     *  - run the resulting executable, unless a compiled binary already exists.
     *
     * @param wilson Name of the Wilson basis.
     * @param model  Model name.
     */
    void compile_run(std::string wilson, std::string model);

    /**
     * @brief Generates the numeric library wrapper for a given (Wilson, model) pair.
     *
     * This step:
     *  - prepares an ::SMParamSetter to compute numeric parameter values,
     *  - creates a ::GeneralNumModelModifier to adapt numeric templates,
     *  - runs a ::NumericTemplateManager / ::CodeGenerator pipeline to
     *    produce the numeric example and parameter files.
     *
     * @param wilson   Wilson basis name.
     * @param model    Model name.
     * @param Q_match  Matching scale (GeV) used later by libraries (currently
     *                 forwarded to compilation stage).
     */
    void generate_numlib(std::string wilson, std::string model);

    /**
     * @brief Generates the numeric wrapper while resolving parameters against
     *        a target model that may differ from the output/cache label.
     *
     * This is needed for SM-like calculations inside a BSM model: generated
     * files keep the ``SM`` label, but the analytic expression may still
     * depend on target-model parameters such as THDM beta.
     */
    void generate_numlib(std::string wilson,
                         std::string output_model,
                         std::string target_model,
                         bool bsm_split_generation = false,
                         bool full_target_generation = false);

    /**
     * @brief Compiles and runs the numeric libraries.
     *
     * Uses ::MakeCompilerStrategy to:
     *  - build the numeric library in the appropriate directory,
     *  - run it using a given matching scale @p Q_match.
     *
     * @param wilson   Wilson basis name.
     * @param model    Model name.
     * @param Q_match  Matching scale passed to the numeric executable.
     */
    void compile_run_libs(std::string wilson, std::string model, double Q_match);

    /**
     * @brief Convenience shortcut to generate, compile, and run the full pipeline.
     *
     * Orchestrates:
     *  - ::generate,
     *  - ::compile_run,
     *  - ::generate_numlib,
     *  - ::compile_run_libs.
     *
     * @param wilson      Wilson basis name.
     * @param model       Model name.
     * @param Q_match     Matching scale for numeric libraries.
     * @param model_path  Path to the C++ model header.
     */
    void calculate(std::string wilson, std::string model, double Q_match, std::string model_path);

    /**
     * @brief Full pipeline with a separate output label and target model.
     */
    void calculate(std::string wilson,
                   std::string output_model,
                   std::string target_model,
                   double Q_match,
                   std::string model_path,
                   bool sm_like_filter,
                  bool bsm_split_generation = false,
                  bool full_target_generation = false);

    /**
     * @brief Retrieves the set of parameter dependencies for a given Wilson basis.
     *
     * The dependencies are discovered via numeric template generation, and
     * correspond to the set of ::InterpretedParam used in the numeric code.
     *
     * @param wilson Wilson basis name.
     * @return Set of interpreted parameters required for this basis.
     *
     * @throws If the @p wilson key is unknown in the internal cache, a
     *         logged error is emitted and the behavior is undefined.
     */
    std::unordered_set<InterpretedParam> get_dependencies(std::string wilson);

    /**
     * @brief Returns the special blocks handled with custom logic.
     *
     * These are the block names for which ::SMParamSetter applies special
     * formulas instead of simply reading values from the parameter proxies.
     *
     * Typical examples include:
     *  - `"KIN"`   : kinetic terms,
     *  - `"WEIN"`  : Weinberg angle,
     *  - `"Finite"`, `"REGPROP"`, `"BETA"`, etc.
     *
     * @return A set containing all special block names.
     */
    std::set<std::string> get_special_blocks();

private:
    /**
     * @brief Checks whether a given binary has already been run/produced.
     *
     * Uses ::stat to test for existence and non-zero size of @p outputBinary.
     *
     * @param outputBinary Path to the binary to test (rvalue to allow moves).
     * @return `true` if the file exists and is non-empty, `false` otherwise.
     */
    bool already_run(std::string&& outputBinary);

    /**
     * @brief Helper to construct the base output binary name.
     *
     * This function is currently unused in the implementation but preserved
     * for convenience.
     *
     * @param wilson Wilson basis name.
     * @param model  Model name.
     * @return A string like `"generated_<wilson>_<model>.cpp"`.
     */
    std::string output_binary_name(std::string& wilson, std::string& model);

    /**
     * @brief Resolve the structural MARTY template index for template models.
     *
     * For the 2HDM this reads MINPAR(24), i.e. the Yukawa type, and uses it to
     * instantiate THDM_Model<N>. Returning std::nullopt means that the model is
     * not a template model.
     */
    std::optional<int> resolve_model_template_index(const std::string& model) const;

    /**
     * @brief Remove stale generated MARTY files when the resolved model signature changes.
     */
    void invalidate_template_model_cache_if_needed(const std::string& wilson,
                                                   const std::string& output_model,
                                                   const std::string& target_model,
                                                   const std::string& model_path,
                                                   std::optional<int> model_template_index,
                                                   bool sm_like_filter,
                                                   bool bsm_split_generation,
                                                   bool full_target_generation) const;

    /// Set of block names requiring special SM handling (e.g. kinematics, angles).
    std::set<std::string> specials_block {"KIN", "WEIN", "Finite", "REGPROP", "BETA"};

    /// Path to the last generated numeric file (reserved for future use).
    std::string num_file_path{};

    /// Cache of dependencies per Wilson basis name.
    std::unordered_map<std::string, std::unordered_set<InterpretedParam>> dependencies;

    /// Core API used to retrieve the current ::Model.
    std::shared_ptr<ICoreAPI<Model>> core_api;

    /// Proxy to SM parameters.
    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> param_proxy_sm;
    
    /// Factory for building interpreter ports (resolvers, mapping adapters).
    std::shared_ptr<IInterpreterPortsFactory> ports;


    /// Optional proxy to BSM parameters (may be null if unused).
    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> param_proxy_bsm = nullptr;
};

#endif
