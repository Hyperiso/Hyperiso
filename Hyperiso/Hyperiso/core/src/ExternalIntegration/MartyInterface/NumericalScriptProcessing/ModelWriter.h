#ifndef MODEL_WRITER_H
#define MODEL_WRITER_H

#include "LineProcessor.h"
#include "ParamWriter.h"
#include <fstream>

class ModelWriter {
private:
    LineProcessor lineProcessor;
    ParamWriter paramwriter;

public:
    ModelWriter(LineProcessor& lineProcessor, ParamWriter& paramwriter);

    void writeModel(std::ifstream& inputFile, std::ofstream& outputFile);

    void writeParam(std::ofstream& paramFile, const std::unordered_map<std::string, double>& params);
};

#endif