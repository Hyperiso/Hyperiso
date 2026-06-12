#include "GeneralModelModifier.h"
#include "MartyRuntimeConfig.h"

GeneralModelModifier::GeneralModelModifier(std::string wilson, std::string model, std::string model_path) {this->wilson = wilson;
        this->model = model;
        std::string model_l = this->model;
        std::transform(model_l.begin(), model_l.end(), model_l.begin(), 
        [](unsigned char c){return std::tolower(c);});
        this->model_path = model_path;

        const auto marty = MartyRuntimeConfig::require_available("GeneralModelModifier");
        if (marty.valid) {
            this->marty_path = marty.marty_header.string();
        }
    }

void GeneralModelModifier::modifyLine(std::string& line) {
    if (line.find("SM_Model sm;") != std::string::npos) {
        bool is_template = ModelFileChecker(this->model_path).isAnyModelTemplate();
        if (is_template) {
            int model_num = 2;
            if(this->model == "THDM"){
                line.replace(line.find("SM_Model"), 8, this->model + "_Model<" + std::to_string(model_num) + ">");
            } else {
                throw std::runtime_error("Unknown model");
            }
        } else {
            line.replace(line.find("SM"), 2, this->model);
        }
    }
    else if (line.find("_SM") != std::string::npos) {
        line.replace(line.find("SM"), 2, this->model);
    }
}

void GeneralModelModifier::addLine(std::ofstream& outputFile, const std::string& currentLine) {
    if (currentLine.find("<iostream>") != std::string::npos) {
        outputFile << currentLine << "\n";
        std::string model_l = this->model;
        std::transform(model_l.begin(), model_l.end(), model_l.begin(),
        [](unsigned char c){return std::tolower(c);});
        outputFile << "#include \"" + this->model_path + "\"" << "\n";
        outputFile << "#include \"" + this->marty_path + "\"" << "\n";
    }
    else {
        outputFile << currentLine << "\n";
    }

}
