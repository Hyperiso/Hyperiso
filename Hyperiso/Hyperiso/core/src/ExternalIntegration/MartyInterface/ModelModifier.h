#ifndef MODEL_MODIFIER_H
#define MODEL_MODIFIER_H

#include <string>
#include <fstream>
#include<iostream>

/**
 * @file ModelModifier.h
 * @brief Declares the base interface for source-level model transformations.
 *
 * This header defines ::ModelModifier, an abstract interface used to
 * modify C++ model templates line by line (e.g. to adapt the model name,
 * inject includes, etc.).
 */

/**
 * @class ModelModifier
 * @ingroup CodeGenerationModule
 * @brief Abstract base class for model source modifiers.
 *
 * A ModelModifier receives each line of a template source file and:
 *  - can mutate the line in-place via ::modifyLine,
 *  - can decide how the line is written via ::addLine.
 *
 * Concrete subclasses (e.g. ::GeneralModelModifier) implement the actual
 * logic to adapt a generic MARTY template to a specific model/Wilson
 * configuration.
 */
class ModelModifier {
public:
    /// Virtual destructor for safe polymorphic use.
    virtual ~ModelModifier() = default;

    /**
     * @brief Modifies a single line of source code in place.
     *
     * Implementations may, for example, replace occurrences of `"SM"`
     * with another model name, or adjust class templates.
     *
     * @param line The line of code to modify.
     */
    virtual void modifyLine(std::string& line) = 0;

    /**
     * @brief Writes a (possibly modified) line to the output file.
     *
     * Default behavior simply writes @p currentLine followed by a newline.
     * Derived classes can override this to inject additional lines before
     * or after the current one, depending on @p addBefore.
     *
     * @param outputFile  Output stream to write to.
     * @param currentLine Line to be written.
     * @param addBefore   Hint that the hook may use to insert content
     *                    before this line (default implementation ignores it).
     */
    virtual void addLine(std::ofstream& outputFile, const std::string& currentLine, bool addBefore) {outputFile << currentLine << "\n";};

protected:
    /// Name of the Wilson operator basis used in the model.
    std::string wilson{};
};

#endif  // MODEL_MODIFIER_H
