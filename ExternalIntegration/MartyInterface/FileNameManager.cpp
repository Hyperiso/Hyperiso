#include "FileNameManager.h"
#include <algorithm>
#include <stdexcept>
#include <iostream>


std::map<std::string, std::shared_ptr<FileNameManager>> FileNameManager::instances;

std::shared_ptr<FileNameManager> FileNameManager::getInstance(const std::string& wilson, const std::string& model) {
    std::string key = wilson + ":" + model;

    auto it = instances.find(key);
    if (it != instances.end()) {
        return it->second;
    }

    if (wilson.empty() || model.empty()) {
        throw std::runtime_error("FileNameManager must be initialized with Wilson and model parameters.");
    }

    std::shared_ptr<FileNameManager> instance(new FileNameManager(wilson, model));
    instances[key] = instance;
    return instance;
}

FileNameManager::FileNameManager(const std::string& wilson, const std::string& model)
    : wilson_(wilson), model_(model) {
    lowercaseWilson_ = toLowercase(wilson_);
    lowercaseModel_ = toLowercase(model_);
    baseDir_ = "libs/" + wilson_ + "_" + model_;
}


std::string FileNameManager::getGeneratedFileName() const {
    return this->root_dir + this->templateDir_ +"generated_" +this->model_ + "_" + this->wilson_ + ".cpp";
}

std::string FileNameManager::getExecutableFileName() const {
    return this->root_dir + this->templateDir_ +"generated_" +this->model_ + "_" + this->wilson_ + "";
}

std::string FileNameManager::getNumGeneratedFileName() const {
    return this->root_dir + this->templateDir_ +"libs/"+this->wilson_+"_"+this->model_+ "/script/example_" + lowercaseWilson_ + "_" + lowercaseModel_ + ".cpp";
}

std::string FileNameManager::getNumExecutableFileName() const {
    return this->root_dir + this->templateDir_ +"libs/"+this->wilson_+"_"+this->model_+ "/bin/example_" + lowercaseWilson_ + "_" + lowercaseModel_ + ".x";
}

std::string FileNameManager::getNumParamFileName() const {
    return this->root_dir + this->templateDir_ +"libs/"+this->wilson_+"_"+this->model_+ "/include/params.h";
}

std::string FileNameManager::getHelperFileName(const std::string &extension) const {
    if(extension == "h")
        return this->root_dir + this->templateDir_ +"libs/"+this->wilson_+"_"+this->model_+ "/include/csv_helper." + extension;
    else if (extension == "cpp")
        return this->root_dir + this->templateDir_ +"libs/"+this->wilson_+"_"+this->model_+ "/src/csv_helper." + extension;
    else
        throw std::runtime_error("Not good, Not good");
}

std::string FileNameManager::getBaseHelperFileName(const std::string &extension) const {
    if(extension == "h")
        return this->root_dir + this->baseTemplateDir_ +"csv_helper." + extension;
    else if (extension == "cpp")
        return this->root_dir + this->baseTemplateDir_ +"csv_helper." + extension;
    else
        throw std::runtime_error("Not good, Not good");
}

std::string FileNameManager::toLowercase(const std::string& str) const {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c) { return std::tolower(c); });
    return result;
}

std::string FileNameManager::getOutputDir() const {
    return this->root_dir + this->templateDir_;
}

std::string FileNameManager::getTemplateDir() const {
    return this->root_dir + this->baseTemplateDir_;
}

std::string FileNameManager::getLibDir() const {
    return this->root_dir + this->templateDir_ + "libs/"+this->wilson_ + "_" + this->model_ + "/";
}