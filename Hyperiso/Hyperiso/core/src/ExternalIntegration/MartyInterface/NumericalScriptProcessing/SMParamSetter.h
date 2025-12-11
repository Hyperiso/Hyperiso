#ifndef SMPARAMSETTER_H
#define SMPARAMSETTER_H

#include <cmath>
#include <set>
#include <string>
#include <iostream>
#include <cmath>

#include "config.hpp"
#include "LhaID.h"
#include "Interpreter.h"
#include "IMartyParameterProxy.h"

/**
 * @file SMParamSetter.h
 * @brief Declares a helper to map Hyperiso parameters into MARTY inputs.
 *
 * SMParamSetter translates internal Hyperiso parameters (blocks, LHA IDs,
 * special cases) into a flat map of name → value suitable for MARTY’s
 * generated code and CSV interfaces.
 */

/**
 * @class SMParamSetter
 * @ingroup CodeGenerationModule
 * @brief Converts Hyperiso parameters into MARTY input parameters.
 *
 * SMParamSetter uses one or two ::IMartyParameterProxy instances to fetch
 * numerical values from SM and BSM parameter sets, applies special rules
 * for certain blocks (e.g. `KIN`, `WEIN`, `REGPROP`, `BETA`), and returns
 * a map of parameter names to numerical values ready for MARTY.
 *
 * Typical usage:
 *  - construct with a model name and set of special blocks,
 *  - repeatedly call ::setParam() for each interpreted parameter.
 */
class SMParamSetter {
public:
    /**
     * @brief Constructs a parameter setter for a given model.
     *
     * @param model          Model name (e.g. `"SM"`, `"THDM"`, `"MSSM"`, `"NMSSM"`).
     * @param special_blocks Set of block names treated with special formulas
     *                       (e.g. `"KIN"`, `"WEIN"`, etc.).
     * @param sm_proxy       Proxy used to access SM parameters.
     * @param bsm_proxy      Optional proxy to access BSM parameters; only used
     *                       if the deduced model type is not ::Model::SM.
     */
    SMParamSetter(const std::string& model, std::set<std::string> special_blocks, std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> sm_proxy, std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> bsm_proxy = nullptr) : special_blocks(special_blocks);

    /**
     * @brief Produces MARTY input values for a given interpreted parameter.
     *
     * The returned map may contain:
     *  - a single entry `name → value` for real parameters,
     *  - two entries `name_rel` and `name_img` for complex parameters,
     *  - special values for certain blocks (handled by ::calculateValue()).
     *
     * @param name             Logical MARTY parameter name.
     * @param interpretedParam Interpreted parameter information (block, LHA ID,
     *                         SM/BSM tag, complex/real).
     * @return Map of generated parameter names to numerical values.
     */
    std::unordered_map<std::string, double> setParam(const std::string& name, const InterpretedParam& interpretedParam);

private:
    /**
     * @brief Computes special values for certain blocks and codes.
     *
     * Handles model-specific logic for blocks such as:
     *  - `"KIN"`         : kinetic term-related quantities,
     *  - `"WEIN"`        : Weinberg angle,
     *  - `"REGPROP"`     : regularization parameters,
     *  - `"BETA"`        : mixing angles in BSM models.
     *
     * @param name             Logical MARTY parameter name.
     * @param interpretedParam Parameter meta-information.
     * @return The computed scalar value.
     */
    scalar_t calculateValue(const std::string& name, const InterpretedParam& interpretedParam);

    Model model_type;   ///< Deduced model type from the string passed to the ctor.

    /// Blocks that should be handled by ::calculateValue() instead of direct lookup.
    std::set<std::string> special_blocks;

    /// Proxy to access SM parameters.
    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> sm_proxy;

    /// Proxy to access BSM parameters when needed.
    std::shared_ptr<IMartyParameterProxy<std::string, LhaID>> bsm_proxy;
};

#endif // SMPARAMSETTER_H
