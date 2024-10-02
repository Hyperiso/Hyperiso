#include "TemplateManager.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "config.hpp"


std::string joinLines(const std::vector<std::string>& lines) {
    std::stringstream ss;
    for (const auto& line : lines) {
        ss << line << "\n";
    }
    return ss.str();
}


void NumericTemplateManager::generateTemplateImpl(const std::string& templateName, const std::string& outputPath) {
    std::string csv_helper= "csv_helper";
    std::string root_path = project_root.data();
    std::string csv_helper_path = root_path + "/DataBase/MartyTemplate/" + csv_helper + ".cpp";
    std::string csv_new_helper_path = templatesDir + "/src/" +csv_helper + ".cpp";
    if (!this->already_generated(csv_new_helper_path)) {
        std::string csv_header_helper_path = root_path + "/DataBase/MartyTemplate/" +csv_helper + ".h";
        std::string csv_header_new_helper_path = templatesDir + "/include/" +csv_helper + ".h";
        std::string command = "cp " + csv_helper_path + " " + csv_new_helper_path;
        std::string command_h = "cp " + csv_header_helper_path + " " + csv_header_new_helper_path;
        system(command.c_str());
        system(command_h.c_str());
    }
    std::string templatePath = templatesDir + "/script/" + templateName + ".cpp";
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


