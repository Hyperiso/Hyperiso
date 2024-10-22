#ifndef FILENAME_MANAGER_H
#define FILENAME_MANAGER_H

#include <string>
#include <memory>
#include <map>
#include "config.hpp"

class FileNameManager {
public:
    static std::shared_ptr<FileNameManager> getInstance(const std::string& wilson = "", const std::string& model = "");

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
    
private:
    FileNameManager(const std::string& wilson, const std::string& model);

    static std::map<std::string, std::shared_ptr<FileNameManager>> instances;

    std::string wilson_;
    std::string model_;
    std::string baseDir_;
    std::string root_dir = project_root.data();
    std::string lowercaseWilson_;
    std::string lowercaseModel_;
    std::string templateDir_ = "/DataBase/MartyWilson/";
    std::string baseTemplateDir_ = "/DataBase/MartyTemplate/";
    std::string toLowercase(const std::string& str) const;
};

#endif // FILENAME_MANAGER_H