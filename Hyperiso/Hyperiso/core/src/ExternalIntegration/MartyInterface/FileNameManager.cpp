#include "FileNameManager.h"

#include <algorithm>
#include <cctype>
#include <filesystem>
#include <stdexcept>
#include <utility>

#include "MartyAdapter.h"
#include "MartyPathProxy.h"

namespace fs = std::filesystem;

void FileNameManager::setTestingRoots(const std::string& templateRoot,
                                      const std::string& baseTemplateRoot,
                                      const std::string& assetsRoot) {
    test_template_root      = templateRoot;
    test_base_template_root = baseTemplateRoot;
    test_assets_root        = assetsRoot;
}

void FileNameManager::clearTestingRoots() {
    test_template_root.clear();
    test_base_template_root.clear();
    test_assets_root.clear();
}

std::shared_ptr<FileNameManager> FileNameManager::getInstance(const std::string& wilson, const std::string& model) {
    if (wilson.empty() || model.empty()) {
        throw std::runtime_error("FileNameManager must be initialized with Wilson and model parameters.");
    }

    auto adapter = std::make_shared<MartyAdapter>();
    auto proxy = std::make_shared<MartyPathProxy>(adapter);
    return std::shared_ptr<FileNameManager>(new FileNameManager(wilson, model, proxy));
}

FileNameManager::FileNameManager(const std::string& wilson,
                                 const std::string& model,
                                 std::shared_ptr<IMartyPathProxy> path_proxy)
    : path_proxy_(std::move(path_proxy)),
      wilson_(wilson),
      model_(model),
      lowercaseWilson_(toLowercase(wilson)),
      lowercaseModel_(toLowercase(model)) {
    if (!path_proxy_) {
        throw std::invalid_argument("FileNameManager requires a non-null IMartyPathProxy.");
    }

    const fs::path output_root = !test_template_root.empty()
        ? fs::path(test_template_root)
        : path_proxy_->get_path(MartyPath::MARTY_TEMP_DIR);

    const fs::path base_template_root = !test_base_template_root.empty()
        ? fs::path(test_base_template_root)
        : path_proxy_->get_path(MartyPath::TEMPLATE_DIR);

    templateDir_ = ensureTrailingSlash(output_root.string());
    baseTemplateDir_ = ensureTrailingSlash(base_template_root.string());

    if (!test_assets_root.empty()) {
        const fs::path mapping_dir = fs::path(test_assets_root) / "input_files" / "marty_mapping";
        smMappingFile_ = (mapping_dir / "sm.json").string();
        bsmMappingFile_ = (mapping_dir / (lowercaseModel_ + ".json")).string();
    } else {
        smMappingFile_ = path_proxy_->get_path(MartyPath::SM_MAPPING_FILE).string();

        const auto bsm_path = path_proxy_->get_optional_path(MartyPath::BSM_MAPPING_FILE);
        if (bsm_path.has_value()) {
            bsmMappingFile_ = bsm_path.value().string();
        }
    }

    fs::create_directories(templateDir_);
}

std::string FileNameManager::toLowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char c){ return static_cast<char>(std::tolower(c)); });
    return value;
}

std::string FileNameManager::ensureTrailingSlash(const std::string& value) {
    if (value.empty()) {
        return value;
    }

    if (value.back() == '/' || value.back() == '\\') {
        return value;
    }

    return value + "/";
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

std::string FileNameManager::getOutputDir() const {
    return templateDir_;
}

std::string FileNameManager::getTemplateDir() const {
    return baseTemplateDir_;
}

std::string FileNameManager::getLibDir() const {
    return templateDir_ + "libs/" + wilson_ + "_" + model_ + "/";
}

std::string FileNameManager::getCsvWilsonFileName() const {
    return templateDir_ + model_ + "_wilson.csv";
}

std::string FileNameManager::getjsondbmodel() const {
    return bsmMappingFile_.value_or(smMappingFile_);
}

std::string FileNameManager::getSmMappingFileName() const {

    return smMappingFile_;
}

std::optional<std::string> FileNameManager::getBsmMappingFileName() const {
    return bsmMappingFile_;
}

std::string FileNameManager::getParamFileName() const {
    return getLibDir() + "bin/paramlist.csv";
}
