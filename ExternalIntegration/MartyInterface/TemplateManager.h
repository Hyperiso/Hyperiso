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

    bool already_generated(const std::string& path) {
        std::ifstream file(path);

        if (!file.is_open()) {
            std::cerr << "Erreur : Impossible d'ouvrir le fichier " << path << std::endl;
            return false;
        }

        std::string firstLine;
        if (std::getline(file, firstLine)) {
            if (firstLine.find("//42") != std::string::npos) {
                return true;
            }
        }

        return false;
    }
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
