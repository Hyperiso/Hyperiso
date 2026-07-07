#ifndef GENERAL_NUM_MODEL_MODIFIER_H
#define GENERAL_NUM_MODEL_MODIFIER_H

#include <map>
#include <string>
#include <unordered_map>

#include "ModelModifier.h"
#include "Extractor.h"
#include "Interpreter.h"
#include "SMParamSetter.h"
#include "ParamWriter.h"
#include "IncludeManager.hpp"
#include "LineProcessor.h"
#include "ModelWriter.h"
#include "FileWriter.h"
#include "FileNameManager.h"

/**
 * @file GeneralNumModelModifier.h
 * @brief Declares the helper used to prepare numeric MARTY models and parameter files.
 *
 * This header defines ::GeneralNumModelModifier, which:
 *  - extracts parameters from a numeric MARTY example,
 *  - interprets them using ::Interpreter and mapping databases,
 *  - computes numerical values via ::SMParamSetter,
 *  - writes both model and parameter files through ::ModelWriter.
 */

/**
 * @class GeneralNumModelModifier
 * @ingroup CodeGenerationModule
 * @brief High-level orchestrator for numeric MARTY model generation.
 *
 * GeneralNumModelModifier coordinates:
 *  - extraction of symbolic parameters from a numeric template
 *    (::Extractor),
 *  - interpretation against mapping databases (::Interpreter),
 *  - conversion into actual SM/BSM numerical values (::SMParamSetter),
 *  - writing of the modified numeric C++ example and parameter CSV
 *    (::ModelWriter).
 *
 * It is dedicated to "numeric" templates that require a parameter file
 * and CSV output for Wilson coefficients.
 */
class GeneralNumModelModifier {
private:
    std::string wilson;     ///< Wilson basis name.
    std::string model;      ///< Model name (e.g. "SM", "THDM", ...).
    bool bsm_split_generation{false}; ///< Numeric wrapper writes extra SM/TOT diagnostics for BSM split libraries.
    bool forceMode; ///< If true, forces rewriting even if marker `//42` is found.
    
    Interpreter                 interpreter;
    std::unique_ptr<SMParamSetter> paramSetter;
    FileWriter                  fileWriter;
    LineProcessor               lineProcessor;
    ModelWriter                 modelWriter;
    ParamWriter                 paramWriter;
    IncludeManager              includeManager;
    
    /// Optional mapping from user names to internal parameter names.
    std::map<std::string, std::string> paramMap;

    /// Resolved numerical values for all parameters.
    std::unordered_map<std::string, double> params;

    /// Map from parameter name to interpreted mapping (block, LhaID, flags).
    std::unordered_map<std::string, InterpretedParam> interpreted_params;


public:
    /**
     * @brief Constructs a numeric model modifier for a given (Wilson, model) pair.
     *
     * @param wilson          Name of the Wilson basis (e.g. "C1").
     * @param model           Name of the model (e.g. "SM", "THDM").
     * @param param_setter    Preconfigured SMParamSetter, moved into this object.
     * @param api             Core API giving access to the active ::Model.
     * @param ports           Factory used to build an ::IParameterResolver for the interpreter.
     * @param force           If true, forces regeneration even when the output file
     *                        appears to be already processed (marker `//42`).
     */
    GeneralNumModelModifier(const std::string& wilson,
                            const std::string& model,
                            std::unique_ptr<SMParamSetter> param_setter,
                            std::shared_ptr<ICoreAPI<Model>> api,
                            std::shared_ptr<IInterpreterPortsFactory> ports,
                            bool force = false)
        : GeneralNumModelModifier(wilson, model, model, std::move(param_setter), std::move(api), std::move(ports), force, false)
    {}

    /**
     * @brief Constructs a numeric model modifier with distinct output and
     *        parameter-resolution models.
     *
     * @param model         Output/cache label used for generated numeric files.
     * @param mapping_model Model whose MARTY mapping should be used to resolve
     *                      symbolic parameters extracted from the generated code.
     */
    GeneralNumModelModifier(const std::string& wilson,
                            const std::string& model,
                            const std::string& mapping_model,
                            std::unique_ptr<SMParamSetter> param_setter,
                            std::shared_ptr<ICoreAPI<Model>> api,
                            std::shared_ptr<IInterpreterPortsFactory> ports,
                            bool force = false,
                            bool bsm_split_generation = false)
        : wilson(wilson),
          model(model),
          bsm_split_generation(bsm_split_generation),
          forceMode(force),
          interpreter(mapping_model, api, ports),
          paramSetter(std::move(param_setter)),
          fileWriter(wilson, model, bsm_split_generation),
          lineProcessor(includeManager, fileWriter, force),
          modelWriter(lineProcessor, paramWriter) {
        initializeParams();
    }

    /**
     * @brief Copy constructor used to duplicate the whole modification pipeline.
     */
    inline GeneralNumModelModifier(const GeneralNumModelModifier& other)
    : wilson(other.wilson),
      model(other.model),
      bsm_split_generation(other.bsm_split_generation),
      forceMode(other.forceMode),
      interpreter(other.interpreter),
      paramSetter(other.paramSetter ? std::make_unique<SMParamSetter>(*other.paramSetter) : nullptr),
      fileWriter(other.wilson, other.model, other.bsm_split_generation),
      lineProcessor(includeManager, fileWriter, other.forceMode),
      modelWriter(lineProcessor, paramWriter),
      paramWriter(other.paramWriter),
      includeManager(other.includeManager),
      paramMap(other.paramMap),
      params(other.params),
      interpreted_params(other.interpreted_params)
{}

    /**
     * @brief Applies the numeric model transformation to a template file.
     *
     * The input numeric example is read from @p inputFile, processed by
     * ::ModelWriter (and ultimately ::LineProcessor), and written to
     * @p outputFile.
     *
     * @param inputFile  Open input stream for the numeric template.
     * @param outputFile Open output stream for the generated C++ file.
     */
    void modify(std::ifstream& inputFile, std::ofstream& outputFile);

    /**
     * @brief Writes the parameter CSV/param file.
     *
     * The given @p params map is passed to ::ModelWriter, which uses
     * ::ParamWriter under the hood.
     *
     * @param paramFile Output stream where the parameter file is written.
     * @param params    Map from parameter name to value.
     */
    void createparamfile(std::ofstream& paramFile, const std::unordered_map<std::string, double> &params);

    /**
     * @brief Returns the interpreted parameter mapping.
     * @return Map from parameter name to ::InterpretedParam.
     */
    std::unordered_map<std::string, InterpretedParam> get_interpreted_param_map() { return this->interpreted_params; }

    /**
     * @brief Returns the map of numerical parameter values.
     * @return Map from parameter name to double value.
     */
    std::unordered_map<std::string, double> get_params() {return this->params;}
private:
    /**
     * @brief Extracts and interprets numeric parameters, then fills the @ref params map.
     *
     * Steps:
     *  - read the param description file via ::FileNameManager::getNumParamFileName,
     *  - extract declarations using ::Extractor,
     *  - interpret them via ::Interpreter,
     *  - convert them into numerical values with ::SMParamSetter and accumulate
     *    into ::params.
     */
    void initializeParams();
};

#endif  // GENERAL_NUM_MODEL_MODIFIER_H