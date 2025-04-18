#pragma once
#include "ParamWriter.h"
#include "IncludeManager.h"
#include "FileWriter.h"

class LineProcessor {
private:
    IncludeManager includeManager;
    FileWriter fileWriter;
    bool done = false;
    bool forceMode;

public:
    LineProcessor(IncludeManager& includeManager, FileWriter& filewriter, bool forceMode)
        : includeManager(includeManager),fileWriter(filewriter), forceMode(forceMode) {}

    void processLine(std::ofstream& outputFile, const std::string& currentLine) {
        if (currentLine.find("//42") != std::string::npos && !forceMode) {
            this->done = true;
        }

        if (done && !forceMode) {
            outputFile << currentLine << "\n";
            return;
        } else if (currentLine.find("int main") != std::string::npos) {
            outputFile << "int main(int argc, char** argv) {\n";
        } 
        else if (currentLine.find("return 0;") != std::string::npos) {
            fileWriter.add_input_reader(outputFile);
            fileWriter.add_argpars(outputFile);
            fileWriter.add_output_writer(outputFile);
            outputFile << currentLine << "\n";
        } else if (currentLine.find("using namespace") != std::string::npos) {
            includeManager.addIncludes(outputFile);
            outputFile << currentLine << "\n";
        } else {
            outputFile << currentLine << "\n";
        }
    }
};