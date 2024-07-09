#include "MartyInterface.h"
#include <memory>

void MartyInterface::run() {
    std::unique_ptr<ModelModifier> smModifier = std::make_unique<SMModelModifier>();
    std::unique_ptr<ModelModifier> thdmModifier = std::make_unique<THDMModelModifier>();

    TemplateManager templateManager("templates");

    templateManager.setModelModifier(std::move(smModifier));
    // templateManager.setModelModifier(std::move(thdmModifier));

    CodeGenerator codeGenerator(templateManager);
    codeGenerator.generate("C7", "generated_main.cpp");

    GppCompilerStrategy compiler;
    compiler.compile("generated_main.cpp", "main");
}
