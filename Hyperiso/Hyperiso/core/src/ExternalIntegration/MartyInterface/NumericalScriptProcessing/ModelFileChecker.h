#ifndef MODEL_FILE_CHECKER_H
#define MODEL_FILE_CHECKER_H

#include <string>
#include <fstream>
#include <iostream>
#include <regex>
#include <optional>
#include <vector>

/**
 * @file ModelFileChecker.h
 * @brief Declares utilities to detect MARTY model classes in C++ files.
 *
 * This header defines ::ModelFileChecker, a helper that inspects a C++
 * source/header file and resolves the C++ class name corresponding to a
 * user-provided MARTY model name.
 */

struct ModelClassInfo {
    std::string class_name;
    bool is_template = false;
};

/**
 * @class ModelFileChecker
 * @ingroup CodeGenerationModule
 * @brief Checks and resolves MARTY model classes in a model header.
 *
 * The resolver intentionally does not force the legacy ``_Model`` suffix.
 * Given a model name such as ``THDM``, it first looks for ``class THDM``, then
 * uppercase and lowercase variants, and only then tries the corresponding
 * ``*_Model`` names.  This keeps existing ``SM_Model`` / ``THDM_Model`` files
 * working while allowing simpler class names in custom MARTY model headers.
 */
class ModelFileChecker {
public:
    /**
     * @brief Constructs a checker for a given file path.
     * @param filePath Path to the C++ source/header file to inspect.
     */
    ModelFileChecker(const std::string& filePath);

    /**
     * @brief Tests whether the file contains any MARTY model template class.
     *
     * This method is preserved for backwards compatibility. New code should
     * prefer ::resolveModelClass.
     */
    bool isAnyModelTemplate() const;

    /**
     * @brief Resolve the concrete C++ class name for a user-provided model name.
     *
     * Candidate order for ``model`` is:
     *   1. ``model``
     *   2. uppercase(``model``)
     *   3. lowercase(``model``)
     *   4. ``model_Model``
     *   5. uppercase(``model``)``_Model``
     *   6. lowercase(``model``)``_Model``
     *
     * Duplicate candidates are ignored. Both non-template and template class
     * declarations are accepted.
     *
     * @throws std::runtime_error if no candidate is found or if the file cannot
     * be opened.
     */
    ModelClassInfo resolveModelClass(const std::string& model) const;

    /**
     * @brief Return the exact candidate sequence used by ::resolveModelClass.
     */
    static std::vector<std::string> modelClassCandidates(const std::string& model);

private:
    std::string readContents() const;
    static bool hasClassDefinition(const std::string& contents,
                                   const std::string& class_name,
                                   bool& is_template);
    static std::string regexEscape(const std::string& value);

    /// Path to the file to inspect.
    std::string filePath;
};

#endif
