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

    std::string paramPath = FileNameManager::getInstance(wilson, model)->getParamFileName();
    std::ofstream paramFile(paramPath);
    numModifier->createparamfile(paramFile);

    if (this->already_generated(outputPath)) {
        return;
    }
    
    std::string tempFilePath = outputPath + ".tmp";
    std::ofstream outputFile(tempFilePath);
    if (!outputFile) {
        return;
    }


    outputFile << "//42" << "\n";
    numModifier->modify(templateFile, outputFile);

    templateFile.close();
    outputFile.close();

    std::string command = "cp " + tempFilePath + " " + outputPath;
    system(command.c_str());
}

void NonNumericTemplateManager::generateTemplateImpl(const std::string& templateName, const std::string& outputPath) {
    std::string templatePath = templatesDir + "/" + templateName + ".cpp";
    std::ifstream templateFile(templatePath);

    if (!templateFile) {
        return;
    }

    if (this->already_generated(outputPath)) {
        return;
    }
    std::ofstream outputFile(outputPath);
    if (!outputFile) {
        return;
    }

    std::string line;
    std::getline(templateFile, line);

    modelModifier->addLine(outputFile, line, true);
    outputFile << "//42" << "\n";

    while (std::getline(templateFile, line)) {
        if (modelModifier) {
            bool addBefore = true;
            modelModifier->modifyLine(line);
            modelModifier->addLine(outputFile, line, addBefore);
        } else {
            outputFile << line << "\n";
        }
    }
}

bool TemplateManagerBase::already_generated(const std::string& path) {
    std::ifstream file(path);

    if (!file.is_open()) {
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
