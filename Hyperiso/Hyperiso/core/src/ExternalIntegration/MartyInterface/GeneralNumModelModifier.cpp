#include "GeneralNumModelModifier.h"

void GeneralNumModelModifier::modify(std::ifstream& inputFile, std::ofstream& outputFile) {
    modelWriter.writeModel(inputFile, outputFile);
}

void GeneralNumModelModifier::createparamfile(std::ofstream& paramFile) {
    modelWriter.writeParam(paramFile);
}

void GeneralNumModelModifier::initializeParams() {
        std::string filename = FileNameManager::getInstance(wilson, model)->getNumParamFileName();
        
        auto extractedParams = Extractor::extract(filename);
        this->interpreted_params = interpreter.interpret(extractedParams);

        for (const auto& [name, interpreted] : this->interpreted_params) {
            paramSetter.setParam(name, interpreted);
        }
    }