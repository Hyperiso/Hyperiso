#ifndef GENERAL_MODEL_MODIFIER_H
#define GENERAL_MODEL_MODIFIER_H

#include <unordered_map>
#include <cctype>
#include <algorithm>
#include <optional>
#include <string>

#include "ModelModifier.h"
#include "ModelFileChecker.h"
#include "config.hpp"

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
    std::string model{};        ///< Target model name.
    std::string model_path{};   ///< Filesystem path to the model header.
    std::string marty_path{};   ///< Path to MARTY's main include.
    std::optional<int> model_template_index{}; ///< Optional template index.
    ModelClassInfo model_class{}; ///< Resolved C++ class and whether it is templated.
    std::string model_instantiation{}; ///< Concrete C++ type written in the generated file.
};

#endif
