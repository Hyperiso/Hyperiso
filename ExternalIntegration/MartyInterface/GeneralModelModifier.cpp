#include "GeneralModelModifier.h"
#include "ModelFileChecker.h"

void GeneralModelModifier::modifyLine(std::string& line) {
    if (line.find("SM_Model sm;") != std::string::npos) {
        line.replace(line.find("SM"), 2, this->model);
        std::cout << line << std::endl;
    }
}

void GeneralModelModifier::addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) {
    if (currentLine.find("marty/models/sm.h") != std::string::npos) {
        outputFile << currentLine << "\n";
        std::string model_l = this->model;
        std::transform(model_l.begin(), model_l.end(), model_l.begin(),
        [](unsigned char c){return std::tolower(c);});
        outputFile << "#include \"../../ExternalIntegration/MARTY/MARTY_INSTALL/include/marty/models/" + model_l +".h\"" << "\n";
    }
    else {
        outputFile << currentLine << "\n";
    }

}