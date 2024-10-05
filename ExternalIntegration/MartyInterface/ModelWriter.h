#pragma once
#include "LineProcessor.h"
#include <fstream>

class ModelWriter {
private:
    LineProcessor lineProcessor;

public:
    ModelWriter(LineProcessor& lineProcessor) : lineProcessor(lineProcessor) {}

    void writeModel(std::ifstream& inputFile, std::ofstream& outputFile) {
        std::string line;
        while (std::getline(inputFile, line)) {
            lineProcessor.processLine(outputFile, line);
        }
    }
};
