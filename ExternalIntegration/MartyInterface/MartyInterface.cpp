#include "MartyInterface.h"
#include <memory>
#include "config.hpp"

void MartyInterface::compile_run(std::string wilson, std::string model) {
    

    GppCompilerStrategy compiler;
    compiler.compile_run("generated_"+wilson+".cpp", "generated_"+wilson);
    this->compiled = true;
}

void MartyInterface::generate(std::string wilson, std::string model) {

    std::unique_ptr<ModelModifier> smModifier;

    if (model == "SM") {
        smModifier = std::make_unique<SMModelModifier>();
    }
    // std::unique_ptr<ModelModifier> thdmModifier = std::make_unique<THDMModelModifier>();

    std::string root_path = project_root.data();
    TemplateManager templateManager(root_path+"/DataBase/MartyTemplate");

    templateManager.setModelModifier(std::move(smModifier));
    // templateManager.setModelModifier(std::move(thdmModifier));

    CodeGenerator codeGenerator(templateManager);
    codeGenerator.generate(wilson, "generated_"+wilson+".cpp");
    this->generated = true;
}
