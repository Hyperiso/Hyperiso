#include "ModelWriter.h"

ModelWriter(LineProcessor& lineProcessor, ParamWriter& paramwriter) : lineProcessor(lineProcessor), paramwriter(paramwriter) {}

void writeModel(std::ifstream& inputFile, std::ofstream& outputFile) {
    std::string line;
    while (std::getline(inputFile, line)) {
        lineProcessor.processLine(outputFile, line);
    }
}

void writeParam(std::ofstream& paramFile, const std::unordered_map<std::string, double>& params) {
    paramwriter.writeParams(paramFile, params);
}