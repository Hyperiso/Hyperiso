#include "MartyInterface.h"
#include <memory>
#include <filesystem>
#include <iostream>
#include "config.hpp"
#include <algorithm>
#include "FileNameManager.h"
#include "GeneralModelModifier.h"

namespace fs = std::filesystem;

std::string to_lowercase(const std::string& str);

void MartyInterface::compile_run(std::string wilson, std::string model) {
    

    GppCompilerStrategy compiler(model, wilson);
    if (!this->already_run(FileNameManager::getInstance(wilson, model)->getNumGeneratedFileName())){
        compiler.compile_run(FileNameManager::getInstance(wilson, model)->getGeneratedFileName(), FileNameManager::getInstance(wilson, model)->getExecutableFileName());
    }
}

void MartyInterface::generate(std::string wilson, std::string model) {

    std::unique_ptr<ModelModifier> smModifier;
    smModifier = std::make_unique<GeneralModelModifier>(wilson, model);

    std::string root_path = project_root.data();
    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NonNumericTemplateManager>(FileNameManager::getInstance(wilson, model)->getTemplateDir());
    templateManager->setModelAndWilson(model, wilson);
    templateManager->setModelModifier(std::move(smModifier));

    CodeGenerator codeGenerator(std::move(templateManager));

    codeGenerator.generate(wilson, FileNameManager::getInstance(wilson, model)->getGeneratedFileName());
}

void MartyInterface::generate_numlib(std::string wilson, std::string model, double Q_match) {
    bool forceMode = false;
    std::unique_ptr<GeneralNumModelModifier> ModelModifier = std::make_unique<GeneralNumModelModifier>(wilson, model, forceMode);
    
    
    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NumericTemplateManager>(FileNameManager::getInstance(wilson, model)->getLibDir());
    templateManager->setModelAndWilson(model, wilson);
    templateManager->setNumModelModifier(std::move(ModelModifier));

    CodeGenerator codeGenerator(std::move(templateManager));

    std::string file_path = FileNameManager::getInstance(wilson, model)->getNumGeneratedFileName();
    codeGenerator.generate(file_path, file_path);
}

void MartyInterface::compile_run_libs(std::string wilson, std::string model, double Q_match) {

    MakeCompilerStrategy compiler(model, wilson);
    compiler.set_Q_match(Q_match);
    compiler.compile_run(FileNameManager::getInstance(wilson, model)->getLibDir(), FileNameManager::getInstance(wilson,model)->getNumExecutableFileName());
}

std::string to_lowercase(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c){return std::tolower(c);});
    return result;
}