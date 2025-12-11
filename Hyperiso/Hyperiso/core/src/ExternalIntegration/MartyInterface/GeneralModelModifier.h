#ifndef GENERAL_MODEL_MODIFIER_H
#define GENERAL_MODEL_MODIFIER_H

#include <unordered_map>
#include <cctype>
#include <algorithm>

#include "ModelModifier.h"
#include "ModelFileChecker.h"
#include "config.hpp"

/**
 * @file GeneralModelModifier.h
 * @brief Declares a general-purpose ModelModifier for non-numeric templates.
 *
 * This class adapts a generic SM-based MARTY template to another model
 * (e.g. THDM) by rewriting class names, model identifiers and injecting
 * the appropriate include directives.
 */

/**
 * @class GeneralModelModifier
 * @ingroup CodeGenerationModule
 * @brief ModelModifier that rewrites SM templates into a target model.
 *
 * Responsibilities:
 *  - Replace occurrences of `SM_Model sm;` with the correct model type,
 *    either as a template (`THDM_Model<2>`) or as a non-template class.
 *  - Replace suffixes like `_SM` with the matching model name.
 *  - Inject the appropriate `#include` lines for the model header and
 *    MARTY's main include.
 */
class GeneralModelModifier : public ModelModifier {
public:
    /**
     * @brief Constructs a modifier for a given (Wilson, model) pair and model header.
     *
     * @param wilson     Name of the Wilson basis (stored in the base class).
     * @param model      Model name (e.g. "THDM").
     * @param model_path Path to the corresponding model header file.
     */
    GeneralModelModifier(std::string wilson, std::string model, std::string model_path);

    /// @copydoc ModelModifier::modifyLine()
    void modifyLine(std::string& line) override;

    /// @copydoc ModelModifier::addLine()
    void addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) override;

private:
    std::string model{};        ///< Target model name.
    std::string model_path{};   ///< Filesystem path to the model header.
    std::string marty_path{};   ///< Path to MARTY's main include.
};

#endif
