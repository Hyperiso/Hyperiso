#ifndef EXTRACTOR_H
#define EXTRACTOR_H

#include <string>
#include <vector>
#include <regex>
#include <fstream>

/**
 * @file Extractor.h
 * @brief Declares utilities to extract parameter declarations from C++ code.
 *
 * This header defines the ::Extractor class, which scans a C++ source file
 * for MARTY-style parameter sanitizer initializations and returns a list of
 * declared parameters.
 */

/**
 * @class Extractor
 * @ingroup CodeGenerationModule
 * @brief Parses C++ code to find MARTY parameter declarations.
 *
 * The Extractor looks for lines of the form:
 * @code
 * csl::InitSanitizer<real_t>    x { "name" };
 * csl::InitSanitizer<complex_t> y { "name" };
 * @endcode
 *
 * and returns a vector of lightweight descriptors containing the variable
 * name, the human-readable name and whether the parameter is complex.
 */
class Extractor {
public:
    /**
     * @struct Parameter
     * @brief Simple descriptor for a discovered parameter.
     *
     * - `type`    : the C++ variable name (e.g. `"x"`),
     * - `name`    : the string name used in the sanitizer (e.g. `"m_b"`),
     * - `complex` : true if declared as `complex_t`, false for `real_t`.
     */
    struct Parameter {
        std::string type;   ///< C++ variable name.
        std::string name;   ///< Logical parameter name from the initializer.
        bool complex;       ///< True for `complex_t`, false for `real_t`.
    };

    /**
     * @brief Extracts all parameter declarations from a given C++ file.
     *
     * The function scans the file with regular expressions looking for
     * `csl::InitSanitizer<real_t>` and `csl::InitSanitizer<complex_t>`
     * patterns.
     *
     * @param filename Path to the C++ file to analyze.
     * @return A vector of discovered ::Parameter descriptors.
     */
    static std::vector<Parameter> extract(const std::string& filename);
};

#endif // EXTRACTOR_H
