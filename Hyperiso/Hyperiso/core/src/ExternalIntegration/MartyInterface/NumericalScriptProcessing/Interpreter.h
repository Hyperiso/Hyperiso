#ifndef INTERPRETER_H
#define INTERPRETER_H

#include <string>
#include <unordered_map>
#include <memory>

#include "Extractor.h"
#include "LhaID.h"
#include "ICoreAPI.h"
#include "IInterpreterPortsFactory.h" 
#include "FileNameManager.h"
#include "InterpretedParam.h"

/**
 * @file Interpreter.h
 * @brief Declares the Interpreter that maps MARTY parameters to Hyperiso IDs.
 *
 * The ::Interpreter reads a set of parameters extracted from MARTY-generated
 * code (via ::Extractor) and, using mapping databases, resolves them to
 * Hyperiso block names and LHA IDs encapsulated in ::InterpretedParam.
 */

/**
 * @class Interpreter
 * @ingroup MappingModule
 * @brief Resolves MARTY parameter names into Hyperiso blocks and IDs.
 *
 * The Interpreter orchestrates:
 *  - the mapping database access (through an ::IInterpreterPortsFactory),
 *  - the current model context (through ::ICoreAPI<Model>),
 *  - the JSON mapping files (located via ::FileNameManager).
 *
 * It takes a list of ::Extractor::Parameter and returns a map of
 * parameter names to ::InterpretedParam, indicating:
 *  - which block they belong to,
 *  - which LHA ID they use,
 *  - whether they are BSM or SM,
 *  - whether they are complex.
 */
class Interpreter {
public:
    /**
     * @brief Constructs an Interpreter for a given model.
     *
     * The constructor builds a parameter resolver using the provided
     * ports factory. It locates the appropriate JSON mapping files
     * for the given model and for the SM using ::FileNameManager.
     *
     * @param model  Model name (e.g. `"SM"`, `"THDM"`, `"MSSM"`).
     * @param api    Core API providing access to the current ::Model value.
     * @param ports  Factory that creates an ::IParameterResolver instance.
     */
    Interpreter(const std::string& model,
                            std::shared_ptr<ICoreAPI<Model>> api,
                            std::shared_ptr<IInterpreterPortsFactory> ports);
    
    /**
     * @brief Interprets a set of extracted parameters.
     *
     * This method:
     *  - queries the underlying resolver,
     *  - forwards the `modelIsSM` flag (depending on ::ICoreAPI<Model>),
     *  - returns an ::InterpretedParam for each input name.
     *
     * @param params Vector of parameters extracted from MARTY code.
     * @return Map of logical parameter name → ::InterpretedParam.
     */
    std::unordered_map<std::string, InterpretedParam>
    interpret(std::vector<Extractor::Parameter>& params);
    
    /// Copy constructor.
    Interpreter(const Interpreter& other);
    /// Copy assignment operator.
    Interpreter& operator=(const Interpreter& other);
    /// Move constructor (defaulted).
    Interpreter(Interpreter&&) noexcept = default;
    /// Move assignment operator (defaulted).
    Interpreter& operator=(Interpreter&&) noexcept = default;

private:
    /// Resolver used to map external names to (block, LHA ID, flags).
    std::unique_ptr<IParameterResolver> resolver;

    /// Accessor to the current model (SM, THDM, SUSY, …).
    std::shared_ptr<ICoreAPI<Model>> marty_api;
};


#endif // INTERPRETER_H