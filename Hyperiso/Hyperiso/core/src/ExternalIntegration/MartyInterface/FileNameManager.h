#ifndef FILENAME_MANAGER_H
#define FILENAME_MANAGER_H

#include <string>
#include <memory>
#include <map>

#include "config.hpp"

/**
 * @file FileNameManager.h
 * @brief Manages file and directory names for generated MARTY code and assets.
 *
 * This header defines ::FileNameManager, a small helper used to build
 * consistent paths for:
 *  - generated C++ sources and executables,
 *  - helper CSV / header files,
 *  - model mapping JSON files,
 *  - library and template directories.
 *
 * The naming is based on a `(wilson, model)` pair and project-wide roots
 * (assets root, template root, etc.).
 */


/**
 * @class FileNameManager
 * @ingroup CodeGenerationModule
 * @brief Centralizes the construction of filenames and directories for generated code.
 *
 * FileNameManager is a small utility class that:
 *  - is instantiated per `(wilson, model)` pair (see ::getInstance()),
 *  - uses build-time variables (e.g. `project_assets_root`) as roots,
 *  - provides a variety of getters for output / template paths.
 *
 * For unit tests, overrides can be provided through ::setTestingRoots().
 */
class FileNameManager {
public:
    /**
     * @brief Retrieves the FileNameManager singleton for a given (wilson, model) pair.
     *
     * Instances are cached in an internal map keyed by `"wilson:model"`.
     * If the instance does not yet exist, it is created.
     *
     * @param wilson  Wilson basis / identifier (must not be empty on first call).
     * @param model   Physics model name (must not be empty on first call).
     * @return Shared pointer to the corresponding FileNameManager instance.
     *
     * @throws std::runtime_error if called with empty @p wilson or @p model
     *         and the instance does not exist yet.
     */
    static std::shared_ptr<FileNameManager> getInstance(const std::string& wilson = "", const std::string& model = "");

    /**
     * @brief Sets custom roots for testing.
     *
     * Allows unit tests to redirect:
     *  - the template directory root,
     *  - the "base" template root,
     *  - the assets root.
     *
     * @param templateRoot      Override for the generated template root.
     * @param baseTemplateRoot  Override for the base (reference) template root.
     * @param assetsRoot        Override for the assets root.
     */
    static void setTestingRoots(const std::string& templateRoot,
                                const std::string& baseTemplateRoot,
                                const std::string& assetsRoot);

    /**
     * @brief Clears any testing roots previously set with ::setTestingRoots().
     */           
    static void clearTestingRoots();

    FileNameManager(const FileNameManager&) = delete;
    FileNameManager& operator=(const FileNameManager&) = delete;
    
    /// @name Generated code and executable names
    ///@{

    /**
     * @brief Returns the path to the main generated C++ file.
     *
     * Pattern: `generated_<model>_<wilson>.cpp` in the template directory.
     */
    std::string getGeneratedFileName() const;

    /**
     * @brief Returns the path to the main generated executable.
     *
     * Pattern: `generated_<model>_<wilson>` in the template directory.
     */
    std::string getExecutableFileName() const;

    /**
     * @brief Returns the path to the numeric example generated C++ file.
     *
     * Pattern:
     *  `libs/<wilson>_<model>/script/example_<wilson>_<model>.cpp`
     *  under the template directory.
     */
    std::string getNumGeneratedFileName() const;

    /**
     * @brief Returns the path to the numeric example executable.
     *
     * Pattern:
     *  `libs/<wilson>_<model>/bin/example_<wilson>_<model>.x`
     *  under the template directory.
     */
    std::string getNumExecutableFileName() const;
    
    ///@}

    /// @name Directories
    ///@{

    /**
     * @brief Returns the output directory where generated files are placed.
     *
     * This is typically the "MartyTemp" directory under the assets root.
     */
    std::string getOutputDir() const;

    /**
     * @brief Returns the base template directory.
     *
     * This is usually the directory that contains the original, unmodified
     * MARTY template files.
     */
    std::string getTemplateDir() const;

    /**
     * @brief Returns the directory for the model-specific library.
     *
     * Pattern:
     *  `<outputDir>/libs/<wilson>_<model>/`
     */
    std::string getLibDir() const;
    
    ///@}

    /// @name Auxiliary / helper files
    ///@{

    /**
     * @brief Returns the path to the numeric parameter header file.
     *
     * Pattern:
     *  `libs/<wilson>_<model>/include/params.h`
     */
    std::string getNumParamFileName() const;

    /**
     * @brief Returns the path to the CSV helper file in the current lib tree.
     *
     * @param extension  File extension, usually `"h"` or `"cpp"`.
     * @return Path to `csv_helper.<extension>` in the generated lib tree.
     *
     * @throws std::runtime_error if the extension is unsupported.
     */
    std::string getHelperFileName(const std::string &extension) const;

    /**
     * @brief Returns the path to the base CSV helper file in the template tree.
     *
     * @param extension  File extension, `"h"` or `"cpp"`.
     * @return Path to `csv_helper.<extension>` under the base template directory.
     *
     * @throws std::runtime_error if the extension is unsupported.
     */
    std::string getBaseHelperFileName(const std::string &extension) const;

    /**
     * @brief Returns the path to the CSV file containing Wilson coefficients.
     *
     * Pattern: `<outputDir>/<model>_wilson.csv`
     */
    std::string getCsvWilsonFileName() const;

    /**
     * @brief Returns the path to the JSON mapping file for the MARTY model.
     *
     * Pattern:
     *  `<assetsRoot>/input_files/marty_mapping/<lowercaseModel>.json`
     */
    std::string getjsondbmodel() const;
    
    /**
     * @brief Returns the path to the parameter list CSV file.
     *
     * Pattern:
     *  `<libDir>/bin/paramlist.csv`
     */
    std::string getParamFileName() const;

    ///@}

private:
    /**
     * @brief Constructs a FileNameManager for a given (wilson, model) pair.
     *
     * Use ::getInstance() instead of calling this constructor directly.
     */
    FileNameManager(const std::string& wilson, const std::string& model);
    
    /// Map of `(wilson:model)` → FileNameManager instances.
    static std::map<std::string, std::shared_ptr<FileNameManager>> instances;
    
    /// Optional testing roots (override normal project-wide paths).
    static inline std::string test_template_root{};
    static inline std::string test_base_template_root{};
    static inline std::string test_assets_root{};

    std::string wilson_;            ///< Original Wilson name.
    std::string model_;             ///< Original model name.
    std::string baseDir_;           ///< Base library directory `libs/<wilson>_<model>`.
    std::string root_dir;           ///< Assets root directory.
    std::string lowercaseWilson_;   ///< Lowercased Wilson name.
    std::string lowercaseModel_;    ///< Lowercased model name.
    std::string templateDir_;       ///< Output template directory.
    std::string baseTemplateDir_;   ///< Base template directory.
    
    /**
     * @brief Converts a string to lowercase.
     * @param str Input string.
     * @return Lowercased copy of @p str.
     */
    std::string toLowercase(const std::string& str) const;
};

#endif
