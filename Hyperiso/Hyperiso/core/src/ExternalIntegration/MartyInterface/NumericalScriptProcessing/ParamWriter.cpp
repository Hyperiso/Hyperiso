#include "ParamWriter.h"

void writeParams(std::ofstream& outputFile, const std::unordered_map<std::string, double>& params) {
    for (const auto& [name, value] : params) {
        std::string real_name = name;
        outputFile << real_name << "," << value << "\n";
    }

}