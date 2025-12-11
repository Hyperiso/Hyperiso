#ifndef MODEL_FILE_CHECKER_H
#define MODEL_FILE_CHECKER_H

#include <string>
#include <fstream>
#include <iostream>
#include <regex>

/**
 * @file ModelFileChecker.h
 * @brief Declares a small utility to detect MARTY model templates in C++ files.
 *
 * This header defines ::ModelFileChecker, a helper that inspects a C++
 * source file and determines whether it contains a MARTY-like template
 * model class definition.
 */

/**
 * @class ModelFileChecker
 * @ingroup CodeGenerationModule
 * @brief Checks if a given C++ file defines a MARTY-style model template.
 *
 * The class scans the file contents and looks for a pattern of the form:
 * `template <...> class Xxx_Model : public Namespace::BaseModel`.
 *
 * It is mostly used to distinguish between generic helper files and
 * actual MARTY model implementations.
 */
class ModelFileChecker {
public:
    /**
     * @brief Constructs a checker for a given file path.
     * @param filePath Path to the C++ source file to inspect.
     */
    ModelFileChecker(const std::string& filePath);

    /**
     * @brief Tests whether the file contains a MARTY model template.
     *
     * This method reads the entire file, removes carriage returns, and
     * runs a regular expression that searches for a `template <...> class
     * *_Model : public Something::Something` pattern.
     *
     * @return `true` if a model template is found, `false` otherwise.
     *
     * @throws std::runtime_error If the file cannot be opened.
     */
    bool isAnyModelTemplate() const;

private:
    /// Path to the file to inspect.
    std::string filePath;
};

#endif