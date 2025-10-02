#include "GeneralNumModelModifier.h"

void GeneralNumModelModifier::modify(std::ifstream& inputFile, std::ofstream& outputFile) {
    modelWriter.writeModel(inputFile, outputFile);
}

void GeneralNumModelModifier::createparamfile(std::ofstream& paramFile, const std::unordered_map<std::string, double> &params) {
    modelWriter.writeParam(paramFile, params);
}

void GeneralNumModelModifier::initializeParams() {
        std::string filename = FileNameManager::getInstance(wilson, model)->getNumParamFileName();
        
        auto extractedParams = Extractor::extract(filename);
        this->interpreted_params = interpreter.interpret(extractedParams);

        for (const auto& [name, interpreted] : this->interpreted_params) {
            auto map_param = paramSetter->setParam(name, interpreted);
            for (const auto& [param_name, param_value] : map_param) {
                params[param_name] = param_value;
            }
        }
    }