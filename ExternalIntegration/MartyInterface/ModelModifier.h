#ifndef MODEL_MODIFIER_H
#define MODEL_MODIFIER_H

#include <string>
#include <fstream>
#include<iostream>
class ModelModifier {
public:
    virtual ~ModelModifier() = default;
    virtual void modifyLine(std::string& line) = 0;
    virtual void addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) {outputFile << currentLine << "\n"; std::cout << "baad" << std::endl;};
};

#endif // MODEL_MODIFIER_H
