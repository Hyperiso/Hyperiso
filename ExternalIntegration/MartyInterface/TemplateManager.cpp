#include "TemplateManager.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

TemplateManager::TemplateManager(const std::string& templatesDir) : templatesDir(templatesDir) {}

void TemplateManager::setModelModifier(std::unique_ptr<ModelModifier> modifier) {
    modelModifier = std::move(modifier);
}

std::string joinLines(const std::vector<std::string>& lines) {
    std::stringstream ss;
    for (const auto& line : lines) {
        ss << line << "\n";
    }
    return ss.str();
}

void TemplateManager::generateTemplate(const std::string& templateName, const std::string& outputPath) {
    std::string templatePath = templatesDir + "/" + templateName + ".cpp";
    std::ifstream templateFile(templatePath);

    if (!templateFile) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier template " << templatePath << std::endl;
        return;
    }

    // Lire tout le contenu du fichier si l'entrée et la sortie sont les mêmes
    std::vector<std::string> lines;
    std::string line;

    std::istringstream stringStream;
    std::istream* input;

    if (templatePath == outputPath) {
        while (std::getline(templateFile, line)) {
            lines.push_back(line);
        }
        templateFile.close();
        stringStream.str(joinLines(lines));
        input = &stringStream;
    } else {
        input = &templateFile;
    }

    std::ofstream outputFile(outputPath);
    if (!outputFile) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier de sortie " << outputPath << std::endl;
        return;
    }

    while (std::getline(*input, line)) {
        if (modelModifier) {
            bool addBefore = true;
            modelModifier->addLine(outputFile, line, addBefore);
            modelModifier->modifyLine(line);
        } else {
            outputFile << line << "\n";
        }
    }
}


