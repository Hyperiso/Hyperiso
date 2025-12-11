#ifndef LINE_PROCESSOR_H
#define LINE_PROCESSOR_H

#include "ParamWriter.h"
#include "IncludeManager.hpp"
#include "FileWriter.h"

/**
 * @file LineProcessor.h
 * @brief Declares a line-based transformer for MARTY-generated C++ files.
 *
 * This header defines the ::LineProcessor class, which reads an existing
 * C++ file line by line and injects additional code fragments (includes,
 * argument parsing, IO calls) at specific markers.
 */

/**
 * @class LineProcessor
 * @ingroup CodeGenerationModule
 * @brief Processes lines of a template C++ file and injects custom code.
 *
 * LineProcessor is typically used to post-process a C++ file generated
 * by MARTY:
 *  - it rewrites the `main` signature to accept `argc, argv`,
 *  - it injects extra `#include` directives,
 *  - it adds argument parsing and parameter IO code before `return 0;`,
 *  - it can stop modifying the file after a special marker (`//42`)
 *    unless `forceMode` is enabled.
 */
class LineProcessor {
private:
    IncludeManager includeManager;  ///< Manager for inserting extra includes.
    FileWriter fileWriter;          ///< Helper to write argument and IO code.
    bool done = false;              ///< Indicates whether further injection is disabled.
    bool forceMode;                 ///< If true, ignore the `//42` marker and keep injecting.

public:
    /**
     * @brief Constructs a LineProcessor with the given helpers and mode.
     *
     * @param includeManager Reference to an ::IncludeManager instance.
     * @param filewriter     Reference to a ::FileWriter instance.
     * @param forceMode      If true, the processor keeps modifying the file
     *                       even after encountering the `//42` marker.
     */
    LineProcessor(IncludeManager& includeManager, FileWriter& filewriter, bool forceMode);

    /**
     * @brief Processes a single line of input and writes to the output file.
     *
     * The behavior depends on the contents of @p currentLine:
     *  - If it contains `"//42"` and `forceMode == false`, further injections
     *    are disabled and subsequent lines are passed through unchanged.
     *  - If it contains `"int main"`, the signature is rewritten as
     *    `int main(int argc, char** argv)`.
     *  - If it contains `"return 0;"`, the method inserts:
     *      - parameter input reader code,
     *      - argument parsing code,
     *      - Wilson output writer code,
     *    before forwarding the `return 0;` line.
     *  - If it contains `"using namespace"`, required includes are injected
     *    and the line is forwarded.
     *  - Otherwise, the line is simply forwarded to @p outputFile.
     *
     * @param outputFile  Stream representing the transformed C++ file.
     * @param currentLine Current line read from the input file.
     */
    void processLine(std::ofstream& outputFile, const std::string& currentLine);
};

#endif