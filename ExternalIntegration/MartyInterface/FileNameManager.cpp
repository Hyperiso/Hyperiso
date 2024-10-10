#include "FileNameManager.h"
#include <algorithm>
#include <stdexcept>

FileNameManager* FileNameManager::instance = nullptr;

FileNameManager* FileNameManager::getInstance(const std::string& wilson, const std::string& model) {
    if (instance == nullptr) {
        if (wilson.empty() || model.empty()) {
            throw std::runtime_error("FileNameManager must be initialized with Wilson and model parameters.");
        }
        instance = new FileNameManager(wilson, model);
    }
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
    return this->root_dir + this->templateDir_ + "script/example_" + lowercaseWilson_ + "_" + lowercaseModel_ + ".cpp";
}

std::string FileNameManager::getNumExecutableFileName() const {
    return this->root_dir + this->templateDir_ + "bin/example_" + lowercaseWilson_ + "_" + lowercaseModel_ + ".x";
}

std::string FileNameManager::getHelperFileName(const std::string &extension) const {
    return this->root_dir + this->templateDir_ + "csv_helper." + extension;
}

std::string FileNameManager::toLowercase(const std::string& str) const {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c) { return std::tolower(c); });
    return result;
}
