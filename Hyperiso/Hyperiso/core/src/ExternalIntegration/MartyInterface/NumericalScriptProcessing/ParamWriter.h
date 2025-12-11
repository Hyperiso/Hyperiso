#ifndef PARAM_WRITER_H
#define PARAM_WRITER_H

#include <unordered_map>
#include <string>
#include <fstream>

#include "config.hpp"
#include "FileNameManager.h"

/**
 * @file ParamWriter.h
 * @brief Declares a small helper to write parameter maps to CSV files.
 *
 * This header defines the ::ParamWriter class, which takes a map of
 * parameter names to numerical values and writes them in a simple
 * `name,value` CSV format.
 */

/**
 * @class ParamWriter
 * @ingroup CodeGenerationModule
 * @brief Helper class to serialize parameter values to a text stream.
 *
 * ParamWriter is typically used in the MARTY code-generation pipeline to dump
 * user parameters or matching inputs into a CSV file that can later be read
 * by a generated executable or an external tool.
 */
class ParamWriter {
public:
    /// Default constructor.
    ParamWriter() =default;

    /**
     * @brief Writes a map of parameters to an output stream.
     *
     * The format of each line is:
     * @code
     * <name>,<value>\n
     * @endcode
     *
     * @param outputFile  Stream to write into (must be opened beforehand).
     * @param params      Map from parameter name to its numeric value.
     */
    void writeParams(std::ofstream& outputFile, const std::unordered_map<std::string, double>& params);

};

#endif