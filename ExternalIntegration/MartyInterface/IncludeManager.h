#pragma once
#include <fstream>

class IncludeManager {
public:
    void addIncludes(std::ofstream& outputFile) {
        outputFile << "#include <fstream>\n";
        outputFile << "#include \"csv_helper.h\"\n";
    }
};