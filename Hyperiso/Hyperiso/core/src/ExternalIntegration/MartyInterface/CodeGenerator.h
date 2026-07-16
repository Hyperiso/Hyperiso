#ifndef CODE_GENERATOR_H
#define CODE_GENERATOR_H

#include <memory>

#include "TemplateManager.h"

/**
 * @file CodeGenerator.h
 * @brief Definition of the CodeGenerator class.
 * 
 * This file defines a class responsible for generating Marty's code using a template manager.
 * The CodeGenerator uses a TemplateManagerBase-derived object to process templates (numerical or analytical part of Marty)
 * and generate output files.
 */

 
/**
 * @class CodeGenerator
 * @brief A class that generates code using a template manager.
 * 
 * The CodeGenerator class utilizes a TemplateManagerBase-derived object (numerical or non-numerical) to generate
 * Marty's code templates based on provided template names and output paths.
 */
class CodeGenerator {
public:
    /**
     * @brief Constructs a CodeGenerator with a unique TemplateManagerBase instance.
     * 
     * @param manager A unique pointer to a TemplateManagerBase-derived object.
     */
    CodeGenerator(std::unique_ptr<TemplateManagerBase> manager);

    /**
     * @brief Generates a code file based on a specified template.
     * 
     * Uses the template manager to process a given Marty's template and generate a Marty file.
     * 
     * @param templateName The name of the template to use.
     * @param outputPath The path where the generated file should be saved.
     */
    void generate(const std::string& templateName, const std::string& outputPath);

private:
    std::unique_ptr<TemplateManagerBase> templateManager; ///< Pointer to the template manager used for code generation.
};

#endif // CODE_GENERATOR_H
