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

void NumericTemplateManager::generateTemplateImpl(const std::string&,
                                                  const std::string& outputPath) {
    namespace fs = std::filesystem;

    auto mgr = FileNameManager::getInstance(this->wilson, this->model);

    const fs::path helper_cpp_src = mgr->getBaseHelperFileName("cpp");
    const fs::path helper_h_src   = mgr->getBaseHelperFileName("h");
    const fs::path helper_cpp_dst = mgr->getHelperFileName("cpp");
    const fs::path helper_h_dst   = mgr->getHelperFileName("h");

    fs::create_directories(helper_cpp_dst.parent_path());
    fs::create_directories(helper_h_dst.parent_path());

    if (!this->already_generated(helper_cpp_dst.string())) {
        fs::copy_file(helper_cpp_src, helper_cpp_dst, fs::copy_options::overwrite_existing);
        fs::copy_file(helper_h_src,   helper_h_dst,   fs::copy_options::overwrite_existing);
    }

    const fs::path templatePath = mgr->getNumGeneratedFileName();
    std::ifstream templateFile(templatePath);

    const fs::path paramPath = mgr->getParamFileName();
    fs::create_directories(paramPath.parent_path());
    std::ofstream paramFile(paramPath);
    std::unordered_map<std::string, double> params = numModifier->get_params();
    numModifier->createparamfile(paramFile, params);

    const fs::path outPathFs = outputPath;
    fs::create_directories(outPathFs.parent_path());

    if (this->already_generated(outputPath)) {
        return;
    }

    const fs::path tmpPath = outPathFs.string() + ".tmp";
    std::ofstream outputFile(tmpPath);
    if (!outputFile) return;

    outputFile << "//42\n";
    numModifier->modify(templateFile, outputFile);
    templateFile.close();
    outputFile.close();

    fs::copy_file(tmpPath, outPathFs, fs::copy_options::overwrite_existing);
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

    modelModifier->addLine(outputFile, line);
    outputFile << "//42" << "\n";

    while (std::getline(templateFile, line)) {
        if (modelModifier) {
            modelModifier->modifyLine(line);
            modelModifier->addLine(outputFile, line);
        } else {
            outputFile << line << "\n";
        }
    }
}

std::unordered_set<InterpretedParam> TemplateManagerBase::get_dependencies() {
    std::unordered_set<InterpretedParam> deps;
    for (auto& [k, v] : numModifier->get_interpreted_param_map()) {
        deps.emplace(v);
    }
    return deps;
}

bool TemplateManagerBase::already_generated(const std::string &path)
{
    std::ifstream file(path);

    if (!file.is_open()) {
        return false;
    }

    // Non-numeric MARTY sources preserve the template's first include and put
    // the generation marker on the following line.  Checking only line one
    // therefore caused every calculate() call to rewrite the analytical source
    // even when its ABI/model/template cache metadata was still valid.
    std::string line;
    for (std::size_t line_number = 0;
         line_number < 16 && std::getline(file, line);
         ++line_number) {
        if (line.find("//42") != std::string::npos) {
            return true;
        }
    }

    return false;
}
