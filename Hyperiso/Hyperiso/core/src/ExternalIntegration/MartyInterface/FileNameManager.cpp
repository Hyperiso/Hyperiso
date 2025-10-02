#include "FileNameManager.h"
#include <algorithm>
#include <stdexcept>

std::map<std::string, std::shared_ptr<FileNameManager>> FileNameManager::instances;

void FileNameManager::setTestingRoots(const std::string& templateRoot,
                                      const std::string& baseTemplateRoot,
                                      const std::string& assetsRoot) {
    test_template_root     = templateRoot;
    test_base_template_root= baseTemplateRoot;
    test_assets_root       = assetsRoot;
}

void FileNameManager::clearTestingRoots() {
    test_template_root.clear();
    test_base_template_root.clear();
    test_assets_root.clear();
}

std::shared_ptr<FileNameManager> FileNameManager::getInstance(const std::string& wilson, const std::string& model) {
    std::string key = wilson + ":" + model;
    if (auto it = instances.find(key); it != instances.end()) return it->second;
    if (wilson.empty() || model.empty()) throw std::runtime_error("FileNameManager must be initialized with Wilson and model parameters.");
    auto inst = std::shared_ptr<FileNameManager>(new FileNameManager(wilson, model));
    instances[key] = inst;
    return inst;
}

FileNameManager::FileNameManager(const std::string& wilson, const std::string& model)
    : wilson_(wilson), model_(model) {

    lowercaseWilson_ = toLowercase(wilson_);
    lowercaseModel_  = toLowercase(model_);
    baseDir_         = "libs/" + wilson_ + "_" + model_;

    const std::string assets = !test_assets_root.empty() ? test_assets_root : std::string(project_assets_root.data());
    const std::string templ  = !test_template_root.empty()? test_template_root : assets + "MartyTemp/";
    const std::string base   = !test_base_template_root.empty()? test_base_template_root : assets + "template/MARTY/";

    root_dir        = assets;
    templateDir_    = templ.back()=='/' ? templ : templ + "/";
    baseTemplateDir_= base.back()=='/' ? base : base + "/";
}

std::string FileNameManager::toLowercase(const std::string& s) const {
    std::string r = s;
    std::transform(r.begin(), r.end(), r.begin(), [](unsigned char c){ return std::tolower(c); });
    return r;
}

std::string FileNameManager::getGeneratedFileName() const {
    return templateDir_ + "generated_" + model_ + "_" + wilson_ + ".cpp";
}
std::string FileNameManager::getExecutableFileName() const {
    return templateDir_ + "generated_" + model_ + "_" + wilson_;
}
std::string FileNameManager::getNumGeneratedFileName() const {
    return templateDir_ + "libs/" + wilson_ + "_" + model_+ "/script/example_" + lowercaseWilson_ + "_" + lowercaseModel_ + ".cpp";
}
std::string FileNameManager::getNumExecutableFileName() const {
    return templateDir_ + "libs/" + wilson_ + "_" + model_+ "/bin/example_" + lowercaseWilson_ + "_" + lowercaseModel_ + ".x";
}
std::string FileNameManager::getNumParamFileName() const {
    return templateDir_ + "libs/" + wilson_ + "_" + model_+ "/include/params.h";
}
std::string FileNameManager::getHelperFileName(const std::string &ext) const {
    if (ext=="h")   return templateDir_ + "libs/" + wilson_ + "_" + model_ + "/include/csv_helper." + ext;
    if (ext=="cpp") return templateDir_ + "libs/" + wilson_ + "_" + model_ + "/src/csv_helper." + ext;
    throw std::runtime_error("Unsupported extension");
}
std::string FileNameManager::getBaseHelperFileName(const std::string &ext) const {
    if (ext=="h" || ext=="cpp") return baseTemplateDir_ + "csv_helper." + ext;
    throw std::runtime_error("Unsupported extension");
}
std::string FileNameManager::getOutputDir() const { return templateDir_; }
std::string FileNameManager::getTemplateDir() const { return baseTemplateDir_; }
std::string FileNameManager::getLibDir() const { return templateDir_ + "libs/" + wilson_ + "_" + model_ + "/"; }
std::string FileNameManager::getCsvWilsonFileName() const { return templateDir_ + model_ + "_wilson.csv"; }
std::string FileNameManager::getjsondbmodel() const { return root_dir + "input_files/marty_mapping/" + lowercaseModel_ + ".json"; }
std::string FileNameManager::getParamFileName() const { return getLibDir() + "bin/paramlist.csv"; }
