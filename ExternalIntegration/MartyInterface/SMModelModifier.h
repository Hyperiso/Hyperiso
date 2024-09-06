#ifndef SM_MODEL_MODIFIER_H
#define SM_MODEL_MODIFIER_H

#include "ModelModifier.h"

class SMModelModifier : public ModelModifier {
public:
    void modifyLine(std::string& line) override;
};

class SMNumModelModifier : public ModelModifier {
public:
    void modifyLine(std::string& line) override {}
    void addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) override;
};
#endif // SM_MODEL_MODIFIER_H
