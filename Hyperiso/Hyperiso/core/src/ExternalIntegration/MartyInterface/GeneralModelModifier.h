#ifndef GENERAL_MODEL_MODIFIER_H
#define GENERAL_MODEL_MODIFIER_H

#include <unordered_map>
#include "ModelModifier.h"
#include "ModelFileChecker.h"
#include "config.hpp"
#include <cctype>
#include <algorithm>

class GeneralModelModifier : public ModelModifier {
public:
    GeneralModelModifier(std::string wilson, std::string model, std::string model_path);
    void modifyLine(std::string& line) override;
    void addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) override;

private:
    std::string model{};
    std::string model_path{};
    std::string marty_path{};
};

#endif // GENERAL_MODEL_MODIFIER_H
