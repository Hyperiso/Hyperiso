#include "TemplateManager.h"

TemplateManager::TemplateManager(const std::string& templatesDir) : templatesDir(templatesDir) {}

void TemplateManager::setModelModifier(std::unique_ptr<ModelModifier> modifier) {
    modelModifier = std::move(modifier);
}

void TemplateManager::generateTemplate(const std::string& templateName, const std::string& outputPath) {
    std::string templatePath = templatesDir + "/" + templateName + ".cpp";
    std::ifstream templateFile(templatePath);
    std::ofstream outputFile(outputPath);
    std::string line;
    while (std::getline(templateFile, line)) {
        if (modelModifier) {
            modelModifier->modifyLine(line);
        }
        outputFile << line << "\n";
    }
}
