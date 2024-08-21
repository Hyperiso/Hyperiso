#include "MartyInterface.h"
#include <memory>
#include "config.hpp"

void MartyInterface::run() {
    

    GppCompilerStrategy compiler;
    compiler.compile("generated_main.cpp", "main");
}

void MartyInterface::generate() {

    std::unique_ptr<ModelModifier> smModifier = std::make_unique<SMModelModifier>();
    // std::unique_ptr<ModelModifier> thdmModifier = std::make_unique<THDMModelModifier>();

    std::string root_path = project_root.data();
    TemplateManager templateManager(root_path+"/DataBase/MartyTemplate");

    templateManager.setModelModifier(std::move(smModifier));
    // templateManager.setModelModifier(std::move(thdmModifier));

    CodeGenerator codeGenerator(templateManager);
    codeGenerator.generate("C7", "generated_main.cpp");
}
