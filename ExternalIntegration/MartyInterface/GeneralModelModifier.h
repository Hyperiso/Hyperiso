#ifndef GENERAL_MODEL_MODIFIER_H
#define GENERAL_MODEL_MODIFIER_H

#include <unordered_map>
#include "ModelModifier.h"
#include "config.hpp"
#include <cctype>
#include <algorithm>

class GeneralModelModifier : public ModelModifier {
public:
    GeneralModelModifier(std::string wilson, std::string model) {this->wilson = wilson;
        this->model = model;
        std::string model_l = this->model;
        std::transform(model_l.begin(), model_l.end(), model_l.begin(), 
        [](unsigned char c){return std::tolower(c);});
        this->model_path = project_root.data() +std::string()+ "/ExternalIntegration/MARTY/MARTY_INSTALL/include/marty/models/" + model_l +".h";
    }
    void modifyLine(std::string& line) override;
    void addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) override;

private:
    std::string model{};
    std::string model_path{};
};

#endif // GENERAL_MODEL_MODIFIER_H
