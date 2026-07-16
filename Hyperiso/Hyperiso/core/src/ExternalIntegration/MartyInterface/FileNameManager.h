#ifndef FILENAME_MANAGER_H
#define FILENAME_MANAGER_H

#include <memory>
#include <optional>
#include <string>

#include "IMartyPathProxy.h"

/**
 * @file FileNameManager.h
 * @brief Manages file and directory names for generated MARTY code and assets.
 *
 * FileNameManager no longer reads build-time project roots directly. It obtains
 * read-only and writable paths through IMartyPathProxy.
 */
class FileNameManager {
public:
    /**
     * @brief Retrieves a FileNameManager for a given (wilson, model) pair.
     *
     * A new instance is returned each time to avoid stale paths after pre-init
     * cache/path overrides.
     */
    static std::shared_ptr<FileNameManager> getInstance(const std::string& wilson = "", const std::string& model = "");

    /**
     * @brief Testing-only root overrides.
     */
    static void setTestingRoots(const std::string& templateRoot,
                                const std::string& baseTemplateRoot,
                                const std::string& assetsRoot);

    static void clearTestingRoots();

    FileNameManager(const FileNameManager&) = delete;
    FileNameManager& operator=(const FileNameManager&) = delete;

    std::string getGeneratedFileName() const;
    std::string getExecutableFileName() const;
    std::string getNumGeneratedFileName() const;
    std::string getNumExecutableFileName() const;

    std::string getOutputDir() const;
    std::string getTemplateDir() const;
    std::string getLibDir() const;

    std::string getNumParamFileName() const;
    std::string getHelperFileName(const std::string &extension) const;
    std::string getBaseHelperFileName(const std::string &extension) const;

    std::string getCsvWilsonFileName() const;

    /**
     * @brief Legacy accessor kept for existing code.
     *
     * Returns the user BSM mapping file when configured, otherwise the read-only
     * SM mapping file. Collision validation is handled by
     * DefaultInterpreterPortsFactory through the existing JSON mapping loader.
     */
    std::string getjsondbmodel() const;

    std::string getSmMappingFileName() const;
    std::optional<std::string> getBsmMappingFileName() const;

    std::string getParamFileName() const;

private:
    FileNameManager(const std::string& wilson,
                    const std::string& model,
                    std::shared_ptr<IMartyPathProxy> path_proxy);

    static std::string toLowercase(std::string value);
    static std::string ensureTrailingSlash(const std::string& value);

    static inline std::string test_template_root{};
    static inline std::string test_base_template_root{};
    static inline std::string test_assets_root{};

    std::shared_ptr<IMartyPathProxy> path_proxy_;

    std::string wilson_;
    std::string model_;
    std::string lowercaseWilson_;
    std::string lowercaseModel_;

    std::string templateDir_;       ///< Writable MARTY temp/output directory.
    std::string baseTemplateDir_;   ///< Read-only MARTY template directory.
    std::string smMappingFile_;
    std::optional<std::string> bsmMappingFile_;
};

#endif
