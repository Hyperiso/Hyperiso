#pragma once
#include "ParamWriter.h"
#include "IncludeManager.h"

class LineProcessor {
private:
    ParamWriter paramWriter;
    IncludeManager includeManager;
    bool done = false;
    bool forceMode;

public:
    LineProcessor(ParamWriter& paramWriter, IncludeManager& includeManager, bool forceMode)
        : paramWriter(paramWriter), includeManager(includeManager), forceMode(forceMode) {}

    void processLine(std::ofstream& outputFile, const std::string& currentLine) {
        if (currentLine.find("//42") != std::string::npos && !forceMode) {
            this->done = true;
        }

        if (done && !forceMode) {
            outputFile << currentLine << "\n";
            return;
        }

        if (currentLine.find("param.") != std::string::npos) {
            std::string paramName = extractParamName(currentLine);
            double paramValue = std::stod(extractParamValue(currentLine));

            if (paramWriter.getParams().find(paramName) != paramWriter.getParams().end()) {
                paramWriter.writeSingleParam(outputFile, paramName, paramValue);
                paramWriter.getParams().erase(paramName);
            } else {
                outputFile << currentLine << "\n";
            }
        } else if (currentLine.find("return 0;") != std::string::npos) {
            paramWriter.writeParams(outputFile);
            outputFile << currentLine << "\n";
        } else if (currentLine.find("using namespace") != std::string::npos) {
            includeManager.addIncludes(outputFile);
            outputFile << currentLine << "\n";
        } else {
            outputFile << currentLine << "\n";
        }
    }

private:
    std::string extractParamName(const std::string& line) {
        size_t startPos = line.find("param.") + 6;
        size_t endPos = line.find(" =");
        return line.substr(startPos, endPos - startPos);
    }

    std::string extractParamValue(const std::string& line) {
        size_t startPos = line.find("= ") + 2;
        size_t endPos = line.find(";");
        return line.substr(startPos, endPos - startPos);
    }
};