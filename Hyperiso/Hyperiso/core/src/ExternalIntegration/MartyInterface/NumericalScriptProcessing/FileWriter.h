
#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <fstream>

#include "FileNameManager.h"

/**
 * @file FileWriter.h
 * @brief Declares utilities to generate C++ code for MARTY wrappers.
 *
 * This header defines the ::FileWriter class, which is responsible for
 * emitting small code fragments into a generated C++ file:
 *  - argument parsing for matching scales,
 *  - input parameter reading,
 *  - Wilson coefficient output writing.
 */

/**
 * @class FileWriter
 * @ingroup CodeGenerationModule
 * @brief Emits C++ helper code to handle CLI, IO and Wilson output.
 *
 * FileWriter writes small C++ snippets into a given `std::ofstream`,
 * tailored for a specific `(wilson, model)` combination. It does not own
 * the stream; the caller is responsible for opening and closing the file.
 */
class FileWriter {
public:
    /**
     * @brief Constructs a FileWriter for a given (wilson, model) pair.
     *
     * @param wilson  Wilson basis or identifier (used in filenames and in code).
     * @param model   Physics model name (used in filenames and in code).
     */
    FileWriter(const std::string& wilson, const std::string& model);

    /**
     * @brief Writes C++ code that handles writing Wilson coefficients to a CSV.
     *
     * The generated code typically:
     *  - computes the output path via ::FileNameManager,
     *  - calls `setMu(Q_match)`,
     *  - calls `writeWilsonCoefficients(...)` with the proper arguments.
     *
     * @param outputFile Stream representing the generated C++ file.
     */
    void add_output_writer(std::ofstream& outputFile);

    /**
     * @brief Writes C++ code to parse command-line arguments.
     *
     * The generated code processes arguments like:
     *  - `--Q_match` / `-Q` to override the matching scale,
     *  - `--help` / `-h` to print a small usage message and exit.
     *
     * @param outputFile Stream representing the generated C++ file.
     */
    void add_argpars(std::ofstream& outputFile);

    /**
     * @brief Writes C++ code to read numerical parameters from disk.
     *
     * The generated code:
     *  - builds the parameter file path via ::FileNameManager,
     *  - opens the file,
     *  - calls `readParams(ParamFile, param.realParams, param.complexParams);`
     *    on a `param_t` object.
     *
     * @param outputFile Stream representing the generated C++ file.
     */
    void add_input_reader(std::ofstream& outputFile);

private:
    std::string model;  ///< Model name associated with the generated code.
    std::string wilson; ///< Wilson basis / identifier.
};

#endif