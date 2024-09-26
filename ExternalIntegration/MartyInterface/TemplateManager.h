#ifndef TEMPLATE_MANAGER_H
#define TEMPLATE_MANAGER_H

#include <string>
#include <fstream>
#include <memory>
#include "ModelModifier.h"

class TemplateManagerBase {
public:
    TemplateManagerBase(const std::string& templatesDir) : templatesDir(templatesDir) {}
    virtual ~TemplateManagerBase() = default;

    virtual void generateTemplate(const std::string& templateName, const std::string& outputPath) = 0;

    void setModelModifier(std::unique_ptr<ModelModifier> modifier) {
        modelModifier = std::move(modifier);
    }

protected:
    std::string templatesDir;
    std::unique_ptr<ModelModifier> modelModifier;
};



class NumericTemplateManager : public TemplateManagerBase {
public:
    NumericTemplateManager(const std::string& templatesDir) : TemplateManagerBase(templatesDir)  {}
    void generateTemplate(const std::string& templateName, const std::string& outputPath) override {
        generateTemplateImpl(templateName, outputPath);
    }

private:
    void generateTemplateImpl(const std::string& templateName, const std::string& outputPath);
};

class NonNumericTemplateManager : public TemplateManagerBase {
public:
    NonNumericTemplateManager(const std::string& templatesDir) : TemplateManagerBase(templatesDir)  {}
    void generateTemplate(const std::string& templateName, const std::string& outputPath) override {
        // Implémentation spécifique pour la gestion des templates non-numériques
        generateTemplateImpl(templateName, outputPath);
    }

private:
    void generateTemplateImpl(const std::string& templateName, const std::string& outputPath);
};


#endif // TEMPLATE_MANAGER_H
