#ifndef CODE_GENERATOR_H
#define CODE_GENERATOR_H

#include "TemplateManager.h"

class CodeGenerator {
public:
    CodeGenerator(std::unique_ptr<TemplateManagerBase> manager);
    void generate(const std::string& templateName, const std::string& outputPath);

private:
    std::unique_ptr<TemplateManagerBase> templateManager;
};

#endif // CODE_GENERATOR_H
