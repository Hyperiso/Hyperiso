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
#include "IncludeManager.hpp"
#include "LineProcessor.h"
#include "ModelWriter.h"
#include "FileWriter.h"
#include "FileNameManager.h"

class GeneralNumModelModifier {
private:
    std::map<std::string, std::string> paramMap;
    std::unordered_map<std::string, double> params;
    std::unordered_map<std::string, InterpretedParam> interpreted_params;
    bool done = false;
    bool forceMode = false;
    int count = 0;
    std::string wilson;
    std::string model;
    Interpreter interpreter;
    std::unique_ptr<SMParamSetter> paramSetter;
    ParamWriter paramWriter;
    FileWriter fileWriter;
    IncludeManager includeManager;
    LineProcessor lineProcessor;
    ModelWriter modelWriter;

public:
    GeneralNumModelModifier(const std::string& wilson, const std::string& model, std::set<std::string>& special_blocks, std::unique_ptr<SMParamSetter> param_setter, std::shared_ptr<ICoreAPI<Model>> api, std::shared_ptr<IInterpreterPortsFactory> ports, bool force = false)
        : wilson(wilson), model(model), forceMode(force), interpreter(model, api, ports), paramSetter(std::move(param_setter)), fileWriter(wilson, model), 
          lineProcessor(includeManager, fileWriter, force), modelWriter(lineProcessor, paramWriter) {
        
        initializeParams();
    }

    inline GeneralNumModelModifier(const GeneralNumModelModifier& other)
    : paramMap(other.paramMap),
      params(other.params),
      interpreted_params(other.interpreted_params),
      done(other.done),
      forceMode(other.forceMode),
      count(other.count),
      wilson(other.wilson),
      model(other.model),
      interpreter(other.interpreter),
      paramSetter(other.paramSetter ? std::make_unique<SMParamSetter>(*other.paramSetter) : nullptr),
      paramWriter(other.paramWriter),
      fileWriter(other.wilson, other.model),
      includeManager(other.includeManager),
      lineProcessor(includeManager, fileWriter, other.forceMode),
      modelWriter(lineProcessor, paramWriter)
{}

    void modify(std::ifstream& inputFile, std::ofstream& outputFile);

    void createparamfile(std::ofstream& paramFile, const std::unordered_map<std::string, double> &params);

    std::unordered_map<std::string, InterpretedParam> get_interpreted_param_map() { return this->interpreted_params; }

    std::unordered_map<std::string, double> get_params() {return this->params;}
private:
    void initializeParams();
};

#endif