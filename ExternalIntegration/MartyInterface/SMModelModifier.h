#ifndef SM_MODEL_MODIFIER_H
#define SM_MODEL_MODIFIER_H

#include <unordered_map>
#include "ModelModifier.h"
#include "Extractor.h"
#include "Interpreter.h"
#include "SMParamSetter.h"

class SMModelModifier : public ModelModifier {
public:
    SMModelModifier(std::string wilson) {this->wilson = wilson;}
    void modifyLine(std::string& line) override;
};

class SMNumModelModifier : public ModelModifier {
public:
    SMNumModelModifier(std::string wilson) : paramSetter(params){this->wilson = wilson;}
    // SMNumModelModifier();
    void modifyLine(std::string& line) override {}
    void addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) override;

private:
    bool done = false;
    int count = 0;
    
    std::unordered_map<std::string, double> params;
    Extractor extractor;
    Interpreter interpreter;
    SMParamSetter paramSetter;

    void processParams(std::ofstream& outputFile);
};
#endif // SM_MODEL_MODIFIER_H
