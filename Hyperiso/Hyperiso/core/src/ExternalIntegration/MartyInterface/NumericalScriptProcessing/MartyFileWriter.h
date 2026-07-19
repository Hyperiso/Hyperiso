
#ifndef MARTY_FILE_WRITER_H
#define MARTY_FILE_WRITER_H

#include <fstream>
#include <string>

#include "FileNameManager.h"

/**
 * @file MartyFileWriter.h
 * @brief Declares utilities to generate C++ code for MARTY wrappers.
 *
 * This header defines the ::MartyFileWriter class, which is responsible for
 * emitting small code fragments into a generated C++ file:
 *  - argument parsing for matching scales,
 *  - input parameter reading,
 *  - Wilson coefficient output writing.
 */

/**
 * @class MartyFileWriter
 * @ingroup CodeGenerationModule
 * @brief Emits C++ helper code to handle CLI, IO and Wilson output.
 *
 * MartyFileWriter writes small C++ snippets into a given `std::ofstream`,
 * tailored for a specific `(wilson, model)` combination. It does not own
 * the stream; the caller is responsible for opening and closing the file.
 */
class MartyFileWriter {
public:
    /**
     * @brief Constructs a MartyFileWriter for a given (wilson, model) pair.
     *
     * @param wilson  Wilson basis or identifier (used in filenames and in code).
     * @param model   Physics model name (used in filenames and in code).
     */
    MartyFileWriter(const std::string& wilson,
                    const std::string& model,
                    bool bsm_split_generation = false,
                    bool full_target_generation = false);

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
    /**
     * @brief Returns true when generated code should call setMu(Q_match).
     *
     * Some Wilson coefficients, starting with C10, must keep the LoopTools
     * dimensional scale at its default value because the matching scale used by
     * HyperIso/SuperIso is not the same object as LoopTools' mudim.
     */
    bool should_set_mudim() const;

    std::string wilson; ///< Wilson basis / identifier.
    std::string model;  ///< Model name associated with the generated code.
    bool bsm_split_generation{false}; ///< Write extra diagnostic columns for BSM-only libraries.
    bool full_target_generation{false}; ///< Treat the split expression as the complete target model.
};

#endif