#pragma once
#include "LineProcessor.h"
#include "ParamWriter.h"
#include <fstream>

class ModelWriter {
private:
    LineProcessor lineProcessor;
    ParamWriter paramwriter;

public:
    ModelWriter(LineProcessor& lineProcessor, ParamWriter& paramwriter) : lineProcessor(lineProcessor), paramwriter(paramwriter) {}

    void writeModel(std::ifstream& inputFile, std::ofstream& outputFile) {
        std::string line;
        while (std::getline(inputFile, line)) {
            lineProcessor.processLine(outputFile, line);
        }
    }

    void writeParam(std::ofstream& paramFile) {
        paramwriter.writeParams(paramFile);
    }
};
