#include "CodeGenerator.h"

CodeGenerator::CodeGenerator(TemplateManager& manager) : templateManager(manager) {}

void CodeGenerator::generate(const std::string& templateName, const std::string& outputPath) {
    templateManager.generateTemplate(templateName, outputPath);
}
