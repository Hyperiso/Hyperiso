#include "GeneralNumModelModifier.h"

void GeneralNumModelModifier::modify(std::ifstream& inputFile, std::ofstream& outputFile) {
    modelWriter.writeModel(inputFile, outputFile);
}

void GeneralNumModelModifier::createparamfile(std::ofstream& paramFile) {
    modelWriter.writeParam(paramFile);
}

void GeneralNumModelModifier::initializeParams() {
        std::string filename = FileNameManager::getInstance(wilson, model)->getNumParamFileName();
        
        auto extractedParams = extractor.extract(filename);
        auto interpretedParams = interpreter.interpret(extractedParams);

        for (const auto& [name, interpreted] : interpretedParams) {
            paramSetter.setParam(name, interpreted);
        }
    }