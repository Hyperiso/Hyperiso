#include "MartyInterface.h"
#include <memory>
#include <filesystem>
#include <iostream>
#include "config.hpp"
#include <algorithm>

namespace fs = std::filesystem;

std::string to_lowercase(const std::string& str);

void MartyInterface::compile_run(std::string wilson, std::string model) {
    

    GppCompilerStrategy compiler;
    if (!this->already_run("libs/" + wilson +"_" + model + "/bin/example_"+ to_lowercase(wilson) +"_"+to_lowercase(model)+".x")){
        compiler.compile_run("generated_"+wilson+".cpp", "generated_"+wilson);
    }
}

void MartyInterface::generate(std::string wilson, std::string model) {

    std::unique_ptr<ModelModifier> smModifier;

    if (model == "SM") {
        smModifier = std::make_unique<SMModelModifier>(wilson);
    }
    // std::unique_ptr<ModelModifier> thdmModifier = std::make_unique<THDMModelModifier>();

    std::string root_path = project_root.data();
    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NonNumericTemplateManager>(root_path + "/DataBase/MartyTemplate");

    templateManager->setModelModifier(std::move(smModifier));
    // templateManager.setModelModifier(std::move(thdmModifier));

    CodeGenerator codeGenerator(std::move(templateManager));
    codeGenerator.generate(wilson, "generated_"+wilson+".cpp");
}

void MartyInterface::generate_numlib(std::string wilson, std::string model) {
    std::unique_ptr<ModelModifier> smModifier;

    if (model == "SM") {
        smModifier = std::make_unique<NumModelModifier>(wilson);
    }
    std::string path = "libs/"+wilson+"_"+model+"/script";
    // std::unique_ptr<ModelModifier> smModifier = std::make_unique<SMNumModelModifier>(wilson);
    std::unique_ptr<TemplateManagerBase> templateManager = std::make_unique<NumericTemplateManager>("libs/"+wilson+"_"+model);

    templateManager->setModelModifier(std::move(smModifier));
    // templateManager.setModelModifier(std::move(thdmModifier));
    fs::path file_path;
    try {

        for (const auto& entry : fs::directory_iterator(path)) {
            if (fs::is_regular_file(entry.path())) {
                file_path = entry.path();
            }
        }

    } catch(const fs::filesystem_error& e) {
        std::cerr << "Erreur: " << e.what() << std::endl;
    }
    CodeGenerator codeGenerator(std::move(templateManager));

    codeGenerator.generate(file_path.stem().string(), file_path.string());
    // this->num_file_path = file_path;
}

void MartyInterface::compile_run_libs(std::string wilson, std::string model) {
    // if(this->num_file_path == "") {
    //     LOG_ERROR("ValueError", "must generate librarie first");
    // }
    MakeCompilerStrategy compiler;
    compiler.compile_run("libs/" + wilson +"_" + model, "/bin/example_"+ to_lowercase(wilson) +"_"+to_lowercase(model)+".x");
}

std::string to_lowercase(const std::string& str) {
    std::string result = str;
    std::transform(result.begin(), result.end(), result.begin(), [](unsigned char c){return std::tolower(c);});
    return result;
}