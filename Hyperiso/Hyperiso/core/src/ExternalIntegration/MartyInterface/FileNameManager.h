// FileNameManager.h
#ifndef FILENAME_MANAGER_H
#define FILENAME_MANAGER_H

#include <string>
#include <memory>
#include <map>
#include "config.hpp"

class FileNameManager {
public:
    static std::shared_ptr<FileNameManager> getInstance(const std::string& wilson = "", const std::string& model = "");

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
    std::string getjsondbmodel() const;

    std::string getParamFileName() const;

private:
    FileNameManager(const std::string& wilson, const std::string& model);

    static std::map<std::string, std::shared_ptr<FileNameManager>> instances;

    static inline std::string test_template_root{};
    static inline std::string test_base_template_root{};
    static inline std::string test_assets_root{};

    std::string wilson_;
    std::string model_;
    std::string baseDir_;
    std::string root_dir;
    std::string lowercaseWilson_;
    std::string lowercaseModel_;
    std::string templateDir_; 
    std::string baseTemplateDir_;
    std::string toLowercase(const std::string& str) const;
};

#endif
