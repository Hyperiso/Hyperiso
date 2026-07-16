#ifndef MODEL_WRITER_H
#define MODEL_WRITER_H

#include <fstream>

#include "LineProcessor.h"
#include "ParamWriter.h"

/**
 * @file ModelWriter.h
 * @brief Declares a helper that writes processed model code and parameters.
 *
 * ::ModelWriter coordinates line-level transformation of a MARTY-generated
 * C++ file (via ::LineProcessor) and exports parameter values using
 * ::ParamWriter into a CSV-like format consumed by the wrapper code.
 */

/**
 * @class ModelWriter
 * @ingroup CodeGenerationModule
 * @brief High-level helper for writing modified model code and parameter files.
 *
 * ModelWriter:
 *  - reads an input C++ model file and passes each line through a
 *    ::LineProcessor to inject IO logic,
 *  - writes parameter values (from a map of name → value) using ::ParamWriter.
 */
class ModelWriter {
private:
    LineProcessor lineProcessor;    ///< Line transformer used for code injection.
    ParamWriter paramwriter;        ///< Writer used to export parameter values.

public:
    /**
     * @brief Constructs a ModelWriter with the given helpers.
     *
     * @param lineProcessor Reference to a configured ::LineProcessor.
     * @param paramwriter   Reference to a ::ParamWriter instance.
     */
    ModelWriter(LineProcessor& lineProcessor, ParamWriter& paramwriter);

    /**
     * @brief Processes the model source and writes the transformed file.
     *
     * Reads @p inputFile line by line and forwards each line through
     * ::LineProcessor::processLine, writing the result to @p outputFile.
     *
     * @param inputFile  Open input stream referring to the original model file.
     * @param outputFile Open output stream where the modified code is written.
     */
    void writeModel(std::ifstream& inputFile, std::ofstream& outputFile);

    /**
     * @brief Writes parameters to the parameter file.
     *
     * Uses ::ParamWriter::writeParams to emit lines of the form
     * `name,value` into @p paramFile.
     *
     * @param paramFile Open output stream for the parameter file.
     * @param params    Map from parameter name to value.
     */
    void writeParam(std::ofstream& paramFile, const std::unordered_map<std::string, double>& params);
};

#endif  // MODEL_WRITER_H