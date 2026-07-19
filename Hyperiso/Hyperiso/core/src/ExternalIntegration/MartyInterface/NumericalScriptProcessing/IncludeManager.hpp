#ifndef INCLUDE_MANAGER_H
#define INCLUDE_MANAGER_H

#include <fstream>

/**
 * @file IncludeManager.h
 * @brief Declares a helper to inject extra includes into generated C++ files.
 *
 * This header defines the ::IncludeManager class, which centralizes the
 * insertion of common `#include` directives (for CSV helpers, IO, etc.)
 * into generated MARTY wrapper code.
 */

/**
 * @class IncludeManager
 * @ingroup CodeGenerationModule
 * @brief Manages insertion of required `#include` directives.
 *
 * IncludeManager writes the standard includes needed by the generated
 * main program that will read parameters and write Wilson coefficients.
 * It does not own the output stream; the caller is responsible for
 * opening and closing the underlying file.
 */
class IncludeManager {
public:
    /**
     * @brief Writes the required `#include` directives to the output stream.
     *
     * Currently this emits:
     *  - `#include <fstream>`
     *  - `#include "csv_helper.h"`
     *
     * @param outputFile The stream representing the C++ file being generated.
     */
    void addIncludes(std::ofstream& outputFile) {
        outputFile << "#include <fstream>\n";
        outputFile << "#include <stdexcept>\n";
        outputFile << "#include <string>\n";
        outputFile << "#include \"csv_helper.h\"\n";
    }
};

#endif