#ifndef TEMPLATE_MANAGER_H
#define TEMPLATE_MANAGER_H

#include <string>
#include <fstream>
#include <memory>
#include "ModelModifier.h"

class TemplateManager {
public:
    TemplateManager(const std::string& templatesDir);
    void setModelModifier(std::unique_ptr<ModelModifier> modifier);
    void generateTemplate(const std::string& templateName, const std::string& outputPath);

private:
    std::string templatesDir;
    std::unique_ptr<ModelModifier> modelModifier;
};

#endif // TEMPLATE_MANAGER_H
