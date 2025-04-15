#include "MartyInterface.h"

namespace fs = std::filesystem;

// std::string to_lowercase(const std::string& str);

void MartyInterface::compile_run(std::string wilson, std::string model) {
    

    GppCompilerStrategy compiler(model, wilson);
    if (!this->already_run(FileNameManager::getInstance(wilson, model)->getNumGeneratedFileName())){
        compiler.compile_run(FileNameManager::getInstance(wilson, model)->getGeneratedFileName(), FileNameManager::getInstance(wilson, model)->getExecutableFileName());
    }
}

void MartyInterface::generate(std::string wilson, std::string model) {

    std::unique_ptr<ModelModifier> smModifier;
    smModifier = std::make_unique<GeneralModelModifier>(wilson, model);

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
    this->dependencies.emplace(wilson, templateManager->get_dependencies());

    CodeGenerator codeGenerator(std::move(templateManager));

    std::string file_path = FileNameManager::getInstance(wilson, model)->getNumGeneratedFileName();
    codeGenerator.generate(file_path, file_path);
}

void MartyInterface::compile_run_libs(std::string wilson, std::string model, double Q_match) {

    MakeCompilerStrategy compiler(model, wilson);
    compiler.set_Q_match(Q_match);
    compiler.compile_run(FileNameManager::getInstance(wilson, model)->getLibDir(), FileNameManager::getInstance(wilson,model)->getNumExecutableFileName());
}

void MartyInterface::calculate(std::string wilson, std::string model, double Q_match, bool new_params) {
    generate(wilson, model);
    compile_run(wilson, model);
    generate_numlib(wilson, model, Q_match);
    compile_run_libs(wilson, model, Q_match);
}

std::unordered_set<Interpreter::InterpretedParam> MartyInterface::get_dependencies(std::string wilson) {
    if (!this->dependencies.contains(wilson)) {
        LOG_ERROR("KeyError", "Trying to access dependencies for unknown wilson coefficient", wilson, "in WilsonInterface.");
    }
    
    return this->dependencies.at(wilson);
}

bool MartyInterface::already_run(std::string&& outputBinary) {
    struct stat buffer;
    if (stat(outputBinary.c_str(), &buffer) != 0) {
        return false;
    }
    if (buffer.st_size == 0) {
        return false;
    }
    std::cout << "Already run !" << std::endl;
    return true;
}

std::string MartyInterface::output_binary_name(std::string& wilson, std::string& model) {
        return "generated_" + wilson+"_" + model + ".cpp";
    }

// std::string to_lowercase(const std::string& str) {
//     std::string result = str;
//     std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c){return std::tolower(c);});
//     return result;
// }