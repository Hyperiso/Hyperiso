#include "TemplateManager.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include "config.hpp"
#include "GeneralNumModelModifier.h"
#include "FileNameManager.h"

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
    std::string csv_helper_path = FileNameManager::getInstance(this->wilson, this->model)->getBaseHelperFileName("cpp");
    std::string csv_new_helper_path = FileNameManager::getInstance(this->wilson, this->model)->getHelperFileName("cpp");
    if (!this->already_generated(csv_new_helper_path)) {
        std::string csv_header_helper_path = FileNameManager::getInstance(this->wilson, this->model)->getBaseHelperFileName("h");
        std::string csv_header_new_helper_path = FileNameManager::getInstance(this->wilson, this->model)->getHelperFileName("h");
        std::string command = "cp " + csv_helper_path + " " + csv_new_helper_path;
        std::string command_h = "cp " + csv_header_helper_path + " " + csv_header_new_helper_path;
        system(command.c_str());
        system(command_h.c_str());
    }
    std::string templatePath = FileNameManager::getInstance(this->wilson, this->model)->getNumGeneratedFileName();
    std::ifstream templateFile(templatePath);

    if (this->already_generated(outputPath)) {
        return;
    }
    
    std::string tempFilePath = outputPath + ".tmp";
    std::ofstream outputFile(tempFilePath);
    if (!outputFile) {
        std::cerr << "Erreur: Impossible d'ouvrir le fichier temporaire " << tempFilePath << std::endl;
        return;
    }

    std::string wilson = "C7";
    bool forceMode = false;
    GeneralNumModelModifier modelModifier(wilson, model, forceMode);

    outputFile << "//42" << "\n";
    std::cout << outputPath << std::endl;
    std::cout << tempFilePath << std::endl;
    modelModifier.modify(templateFile, outputFile);

    templateFile.close();
    outputFile.close();

    std::string command = "cp " + tempFilePath + " " + outputPath;
    system(command.c_str());
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
            modelModifier->modifyLine(line);
        } else {
            outputFile << line << "\n";
        }
    }
}