#ifndef TEMPLATE_MANAGER_H
#define TEMPLATE_MANAGER_H

#include <string>
#include <fstream>
#include <memory>
#include "ModelModifier.h"
#include "GeneralNumModelModifier.h"

class TemplateManagerBase {
public:
    TemplateManagerBase(const std::string& templatesDir) : templatesDir(templatesDir) {}
    virtual ~TemplateManagerBase() = default;

    void generateTemplate(const std::string& templateName, const std::string& outputPath) {
        generateTemplateImpl(templateName, outputPath);
    }

    void setModelModifier(std::unique_ptr<ModelModifier> modifier) { modelModifier = std::move(modifier); }

    void setNumModelModifier(std::unique_ptr<GeneralNumModelModifier> modifier) { numModifier = std::move(modifier); }

    void setModelAndWilson(std::string model, std::string wilson) { this->model = model; this->wilson = wilson; }

    std::unordered_set<InterpretedParam> get_dependencies();

protected:
    std::string templatesDir;
    std::string wilson;
    std::string model;
    std::unique_ptr<ModelModifier> modelModifier;
    std::unique_ptr<GeneralNumModelModifier> numModifier;

    virtual void generateTemplateImpl(const std::string& templateName, const std::string& outputPath) = 0;
    bool already_generated(const std::string& path);
};



class NumericTemplateManager : public TemplateManagerBase {
public:
    NumericTemplateManager(const std::string& templatesDir) : TemplateManagerBase(templatesDir)  {}

private:
    void generateTemplateImpl(const std::string& templateName, const std::string& outputPath) override;
};

class NonNumericTemplateManager : public TemplateManagerBase {
public:
    NonNumericTemplateManager(const std::string& templatesDir) : TemplateManagerBase(templatesDir)  {}

private:
    void generateTemplateImpl(const std::string& templateName, const std::string& outputPath) override;
};


#endif
