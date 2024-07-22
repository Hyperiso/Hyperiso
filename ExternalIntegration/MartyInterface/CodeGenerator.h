#ifndef CODE_GENERATOR_H
#define CODE_GENERATOR_H

#include "TemplateManager.h"

class CodeGenerator {
public:
    CodeGenerator(TemplateManager& manager);
    void generate(const std::string& templateName, const std::string& outputPath);

private:
    TemplateManager& templateManager;
};

#endif // CODE_GENERATOR_H
