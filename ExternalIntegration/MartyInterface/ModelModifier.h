#ifndef MODEL_MODIFIER_H
#define MODEL_MODIFIER_H

#include <string>

class ModelModifier {
public:
    virtual ~ModelModifier() = default;
    virtual void modifyLine(std::string& line) = 0;
    virtual void addLine(std::ofstream& outputFile, const std::string& currentLine) = 0;
};

#endif // MODEL_MODIFIER_H
