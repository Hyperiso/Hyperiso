#include "TemplateManager.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>



std::string joinLines(const std::vector<std::string>& lines) {
    std::stringstream ss;
    for (const auto& line : lines) {
        ss << line << "\n";
    }
    return ss.str();
}


void NumericTemplateManager::generateTemplateImpl(const std::string& templateName, const std::string& outputPath) {
    std::string templatePath = templatesDir + "/" + templateName + ".cpp";
    std::ifstream templateFile(templatePath);
    std::vector<std::string> lines;
    std::string line;
    if (!templateFile) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier template " << templatePath << std::endl;
        return;
    }

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
    
    if (this->already_generated(outputPath)) {
        return;
    }
    
    std::ofstream outputFile(outputPath);
    if (!outputFile) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier de sortie " << outputPath << std::endl;
        return;
    }

    outputFile << "//42" << "\n"; 
    while (std::getline(*input, line)) {
        if (modelModifier) {
            bool addBefore = true;
            modelModifier->addLine(outputFile, line, addBefore);
            modelModifier->modifyLine(line);  // Modification de la ligne si nécessaire
        } else {
            outputFile << line << "\n";
        }
    }
}

void NonNumericTemplateManager::generateTemplateImpl(const std::string& templateName, const std::string& outputPath) {
    std::string templatePath = templatesDir + "/" + templateName + ".cpp";
    std::ifstream templateFile(templatePath);

    if (!templateFile) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier template " << templatePath << std::endl;
        return;
    }

    if (this->already_generated(outputPath)) {
        return;
    }
    std::ofstream outputFile(outputPath);
    if (!outputFile) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier de sortie " << outputPath << std::endl;
        return;
    }

    std::string line;
    std::getline(templateFile, line);

    modelModifier->addLine(outputFile, line, true);
    outputFile << "//42" << "\n";

    while (std::getline(templateFile, line)) {
        if (modelModifier) {
            bool addBefore = true;
            modelModifier->addLine(outputFile, line, addBefore);
            modelModifier->modifyLine(line);  // Modification de la ligne si nécessaire
        } else {
            outputFile << line << "\n";
        }
    }
}

// void TemplateManager::generateTemplate(const std::string& templateName, const std::string& outputPath) {
//     std::string templatePath = templatesDir + "/" + templateName + ".cpp";
//     std::ifstream templateFile(templatePath);

//     if (!templateFile) {
//         std::cerr << "Erreur: Impossible d'ouvrir le fichier template " << templatePath << std::endl;
//         return;
//     }

//     // Lire tout le contenu du fichier si l'entrée et la sortie sont les mêmes
//     std::vector<std::string> lines;
//     std::string line;

//     std::istringstream stringStream;
//     std::istream* input;

//     if (templatePath == outputPath) {
//         while (std::getline(templateFile, line)) {
//             lines.push_back(line);
//         }
//         templateFile.close();
//         stringStream.str(joinLines(lines));
//         input = &stringStream;
//     } else {
//         input = &templateFile;
//     }

//     std::ofstream outputFile(outputPath);
//     if (!outputFile) {
//         std::cerr << "Erreur: Impossible d'ouvrir le fichier de sortie " << outputPath << std::endl;
//         return;
//     }

//     while (std::getline(*input, line)) {
//         if (modelModifier) {
//             bool addBefore = true;
//             modelModifier->addLine(outputFile, line, addBefore);
//             modelModifier->modifyLine(line);
//         } else {
//             outputFile << line << "\n";
//         }
//     }
// }


