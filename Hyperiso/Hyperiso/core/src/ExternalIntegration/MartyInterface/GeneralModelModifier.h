#ifndef GENERAL_MODEL_MODIFIER_H
#define GENERAL_MODEL_MODIFIER_H

#include <unordered_map>
#include <cctype>
#include <algorithm>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "ModelModifier.h"
#include "ModelFileChecker.h"
#include "config.hpp"
#include "GeneralEnum.h"

/**
 * @file GeneralModelModifier.h
 * @brief Declares a general-purpose ModelModifier for non-numeric templates.
 *
 * This class adapts a generic SM-based MARTY template to another model
 * by resolving the model class declared in the user-provided model file and
 * rewriting the default ``SM_Model sm;`` instantiation accordingly.
 */

/**
 * @class GeneralModelModifier
 * @ingroup CodeGenerationModule
 * @brief ModelModifier that rewrites SM templates into a target model.
 *
 * Responsibilities:
 *  - Replace occurrences of ``SM_Model sm;`` with the resolved target model
 *    class, optionally instantiated with a template index when the target class
 *    is templated.
 *  - Replace suffixes like ``_SM`` with the matching model name.
 *  - Inject the appropriate ``#include`` lines for the model header and
 *    MARTY's main include.
 *
 * The model class resolver does not force the historical ``_Model`` suffix.
 * Given ``mty_model_name=THDM``, it first looks for ``class THDM``, then case
 * variants, and only then tries ``THDM_Model`` and its case variants.
 *
 * For SM-like generation inside a BSM model, the modifier can keep the output
 * suffix ``_SM`` while instantiating the target model class and inserting a
 * MARTY filter that disables every particle not present in ``SM_Model``.
 */
class GeneralModelModifier : public ModelModifier {
public:
    /**
     * @brief Constructs a modifier for a given (Wilson, model) pair and model header.
     *
     * @param wilson     Name of the Wilson basis (stored in the base class).
     * @param model      Model name (e.g. "THDM").
     * @param model_path Path to the corresponding model header file.
     * @param model_template_index Optional template index used when the resolved
     *        model class is templated, e.g. ``THDM_Model<2>``.
     */
    GeneralModelModifier(std::string wilson,
                         std::string model,
                         std::string model_path,
                         std::optional<int> model_template_index = std::nullopt);

    /**
     * @brief Constructs a modifier with a separate output label and target model.
     *
     * @param wilson     Name of the Wilson basis.
     * @param output_model Model label used in generated libraries/files, e.g. ``SM``.
     * @param target_model Model class name to instantiate, e.g. ``THDM``.
     * @param model_path Path to the target model header.
     * @param model_template_index Optional template index for target templated models.
     * @param disable_non_sm_particles If true, insert a MARTY filter disabling all
     *        particles of the target model that are not present in ``SM_Model``.
     */
    GeneralModelModifier(std::string wilson,
                         std::string output_model,
                         std::string target_model,
                         std::string model_path,
                         std::optional<int> model_template_index,
                         bool disable_non_sm_particles,
                         bool bsm_split_generation = false,
                         bool full_target_generation = false,
                         bool tree_first_fallback = false,
                         MartyOrderPolicy order_policy = MartyOrderPolicy::AUTO);

    /// @copydoc ModelModifier::modifyLine()
    void modifyLine(std::string& line) override;

    /// @copydoc ModelModifier::addLine()
    void addLine(std::ofstream& outputFile, const std::string& currentLine) override;

    /**
     * @brief Resolve the C++ model instantiation used by generated MARTY code.
     *
     * Examples: ``SM``, ``SM_Model``, ``THDM<2>``, ``THDM_Model<2>``.
     */
    static std::string resolveModelInstantiation(const std::string& model,
                                                 const std::string& model_path,
                                                 std::optional<int> model_template_index = std::nullopt);

    /**
     * @brief Build the cache signature written in generated analytical files.
     */
    static std::string modelSignature(const std::string& model,
                                      const std::string& model_path,
                                      std::optional<int> model_template_index = std::nullopt);

private:
    static std::string makeSmFilterHelper();
    static std::string makeTreeLevelWilsonHelper();
    static void replaceWilsonCallWithHelper(std::string& line);
    std::string orderPolicyPreamble() const;
    bool usesRegPropSplit() const;
    bool usesGenericTreeFirst() const;
    static void replaceWilsonOrderArgument(std::string& line);
    bool consumeTreeSafeWilsonCall(std::ofstream& outputFile,
                                   const std::string& currentLine,
                                   bool pair_return,
                                   bool count_graphs);
    void emitTreeSafeWilsonCall(std::ofstream& outputFile,
                                bool pair_return,
                                bool count_graphs);

    std::string output_model{}; ///< Model label used in generated file/library names.
    std::string target_model{}; ///< Model class name to instantiate.
    std::string model_path{};   ///< Filesystem path to the model header.
    std::string marty_path{};   ///< Path to MARTY's main include.
    std::optional<int> model_template_index{}; ///< Optional template index.
    bool disable_non_sm_particles{false}; ///< Whether to add the SM-like filter.
    bool bsm_split_generation{false}; ///< Use the dedicated split-reg_prop generation or BSM-only filter.
    bool full_target_generation{false}; ///< Keep the complete target-model expression instead of filtering to BSM diagrams.
    bool tree_first_fallback{false}; ///< For loop-only templates, test TreeLevel before evaluating OneLoop.
    MartyOrderPolicy order_policy{MartyOrderPolicy::AUTO}; ///< Explicit BSM MARTY order policy.
    bool inside_calculate_function{false}; ///< Internal line-rewrite state for BSM split mode.
    bool skip_old_main{false}; ///< Internal line-rewrite state for BSM split mode.
    bool expression_returned{false}; ///< Whether the calculation body already returned its primary expression.
    std::string generic_builder_name{}; ///< Rewritten calculate function name for generic tree-first mode.
    bool buffering_tree_safe_wilson_call{false}; ///< Buffer a multi-line Wilson call before emitting the safe tree probe.
    std::vector<std::string> tree_safe_wilson_call_lines{}; ///< Buffered computeWilsonCoefficients call.
    ModelClassInfo model_class{}; ///< Resolved C++ class and whether it is templated.
    std::string model_instantiation{}; ///< Concrete C++ type written in the generated file.
};

#endif
