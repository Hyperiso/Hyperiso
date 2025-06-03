#include "GeneralModelModifier.h"

GeneralModelModifier::GeneralModelModifier(std::string wilson, std::string model, std::string model_path) {this->wilson = wilson;
        this->model = model;
        std::string model_l = this->model;
        std::transform(model_l.begin(), model_l.end(), model_l.begin(), 
        [](unsigned char c){return std::tolower(c);});
        this->model_path = model_path;
        // this->model_path = project_tp_root.data() + std::string() + "MARTY/MARTY_INSTALL/include/marty/models/" + model_l +".h";
        this->marty_path = project_tp_root.data() + std::string() + "MARTY/MARTY_INSTALL/include/marty.h";
    }

void GeneralModelModifier::modifyLine(std::string& line) {
    if (line.find("SM_Model sm;") != std::string::npos) {
        bool is_template = ModelFileChecker(this->model_path).isAnyModelTemplate();
        if (is_template) {
            int model_num = 1;
            // std::cout << this->model << std::endl;
            if(this->model == "THDM"){
                line.replace(line.find("SM_Model"), 8, this->model + "_Model<" + std::to_string(model_num) + ">");
            } else {
                std::runtime_error("Unknown model");
            }
        } else {
            line.replace(line.find("SM"), 2, this->model);
        }
    }
    else if (line.find("_SM") != std::string::npos) {
        line.replace(line.find("SM"), 2, this->model);
    }
}

void GeneralModelModifier::addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) {
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