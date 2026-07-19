#include "ParamWriter.h"

#include <iomanip>

void ParamWriter::writeParams(std::ofstream& outputFile, const std::unordered_map<std::string, double>& params) {
    // Keep the historical MARTY parameter-file representation.  C9 in the
    // THDM is a small cancellation-sensitive quantity, and changing the CSV
    // serialization from the legacy six significant digits to full binary
    // precision changes the established numerical result even though the
    // analytical Wilson expression is identical.
    outputFile << std::defaultfloat << std::setprecision(6);

    for (const auto& [name, value] : params) {
        std::string real_name = name;
        outputFile << real_name << "," << value << "\n";
    }

}