// #include "ModelModifier.h"
#include <map>
#include <string>
#include <unordered_map>
#include "Extractor.h"
#include "Interpreter.h"
#include "SMParamSetter.h"
#include "ParamWriter.h"
#include "IncludeManager.h"
#include "LineProcessor.h"
#include "ModelWriter.h"
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
    IncludeManager includeManager;
    LineProcessor lineProcessor;
    ModelWriter modelWriter;

public:
    GeneralNumModelModifier(const std::string& wilson, const std::string& model, bool force = false)
        : wilson(wilson), model(model), forceMode(force), paramSetter(params), paramWriter(params, wilson), 
          lineProcessor(paramWriter, includeManager, force), modelWriter(lineProcessor) {
        
        initializeParams();
    }

    void modify(std::ifstream& inputFile, std::ofstream& outputFile) {
        modelWriter.writeModel(inputFile, outputFile);
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
