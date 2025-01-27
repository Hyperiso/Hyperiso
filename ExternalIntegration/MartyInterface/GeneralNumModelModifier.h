#ifndef GENERAL_NUM_MODEL_MODIFIER_H
#define GENERAL_NUM_MODEL_MODIFIER_H

#include <map>
#include <string>
#include <unordered_map>
#include "ModelModifier.h"
#include "Extractor.h"
#include "Interpreter.h"
#include "SMParamSetter.h"
#include "ParamWriter.h"
#include "IncludeManager.h"
#include "LineProcessor.h"
#include "ModelWriter.h"
#include "FileWriter.h"
#include "FileNameManager.h"

class GeneralNumModelModifier {
private:
    std::map<std::string, std::string> paramMap;
    std::unordered_map<std::string, double> params;
    bool done = false;
    bool forceMode = false;
    int count = 0;
    std::string wilson;
    std::string model;
    Extractor extractor;
    Interpreter interpreter;
    SMParamSetter paramSetter;
    ParamWriter paramWriter;
    FileWriter fileWriter;
    IncludeManager includeManager;
    LineProcessor lineProcessor;
    ModelWriter modelWriter;

public:
    GeneralNumModelModifier(const std::string& wilson, const std::string& model, bool force = false)
        : wilson(wilson), model(model), forceMode(force), interpreter(model), paramSetter(params, model), fileWriter(wilson, model), paramWriter(params, wilson, model), 
          lineProcessor(includeManager, fileWriter, force), modelWriter(lineProcessor, paramWriter) {
        
        initializeParams();
    }

    void modify(std::ifstream& inputFile, std::ofstream& outputFile) {
        modelWriter.writeModel(inputFile, outputFile);
    }

    void createparamfile(std::ofstream& paramFile) {
        modelWriter.writeParam(paramFile);
    }

private:
    void initializeParams() {
        std::string filename = FileNameManager::getInstance(wilson, model)->getNumParamFileName();
        
        auto extractedParams = extractor.extract(filename);
        auto interpretedParams = interpreter.interpret(extractedParams);

        for (const auto& [name, interpreted] : interpretedParams) {
            paramSetter.setParam(name, interpreted);
        }
    }
};

#endif