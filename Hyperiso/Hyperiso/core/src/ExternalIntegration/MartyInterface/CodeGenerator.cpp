#include "CodeGenerator.h"

CodeGenerator::CodeGenerator(std::unique_ptr<TemplateManagerBase> manager) : templateManager(std::move(manager)) {}

void CodeGenerator::generate(const std::string& templateName, const std::string& outputPath) {
    templateManager->generateTemplate(templateName, outputPath);
}
