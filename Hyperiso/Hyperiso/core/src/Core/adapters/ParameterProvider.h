#ifndef PARAMETERPROVIDER_H
#define PARAMETERPROVIDER_H

#include "IDataProvider.h"
#include "Include.h"
#include "Parameters.h"

/**
 * @file ParameterProvider.h
 * @brief High-level access to parameter values and uncertainties.
 *
 * This header declares the ParameterProvider class, a concrete implementation
 * of the generic IDataProvider interface. It provides a unified, callable way
 * to retrieve:
 *  - central values of parameters,
 *  - statistical uncertainties,
 *  - systematic uncertainties,
 *  - combined uncertainties,
 * as well as existence checks and scales for parameter blocks.
 *
 * ParameterProvider is built on top of the Parameters singleton infrastructure
 * and uses ParamId / (block, LhaID) pairs to identify parameters.
 */

/**
 * @class ParameterProvider
 * @ingroup DataProvidersModule
 * @brief Provides access to parameter values, errors, and existence checks.
 *
 * ParameterProvider wraps the Parameters subsystem and exposes a convenient
 * callable interface:
 *
 *  - `operator()(ParamId, DataType)` to fetch central value or uncertainties
 *    of a parameter identified by ParamId.
 *  - `operator()(block, LhaID, DataType)` to fetch the same information using
 *    the usual SLHA-style (block, index) access.
 *
 * It can optionally be “bound” to a given ParameterType (e.g. SM, BSM) at
 * construction time. In that case, calls using `(block, LhaID, ...)` will
 * implicitly use that parameter type.
 *
 * Typical usage:
 * @code
 * ParameterProvider p_sm(ParameterType::SM);
 * ParamId pid{ParameterType::SM, "MASS", LhaID(25)};
 *
 * double mh      = p_sm(pid);                                     // VALUE
 * double sigma_s = p_sm(pid, DataType::STD_STAT);                 // stat. error
 * double sigma_t = p_sm("MASS", LhaID(25), DataType::STD_COMBINED);
 * bool   exists  = p_sm.exists(pid);
 * @endcode
 */
class ParameterProvider : public IDataProvider<ParameterProvider> {
public:
    /**
     * @brief Default constructor (untyped provider).
     *
     * When constructed without a ParameterType, the provider can still be used
     * with ParamId-based access. Block-based access `(block, LhaID, ...)`
     * requires an explicit type (see the other constructor).
     */
    ParameterProvider() = default;

    /**
     * @brief Constructs a provider bound to a specific ParameterType.
     *
     * Ensures that the corresponding Parameters instance exists, and stores
     * the type internally for later block-based queries.
     *
     * @param p_type ParameterType handled by this provider (e.g. SM, BSM).
     */
    inline ParameterProvider(ParameterType p_type) : p_type(p_type) { Parameters::GetInstance(p_type); }

    /**
     * @brief Retrieves a parameter value based on ParamId and requested data type.
     *
     * The ParamId fully encodes the parameter type, block, and LHA identifier.
     *
     * @param pid    The parameter ID.
     * @param d_type Type of data requested (central value, stat error, syst error, etc.).
     * @return The requested scalar value.
     */
    scalar_t operator()(const ParamId& pid, DataType d_type=DataType::VALUE) const;

    /**
     * @brief Retrieves a parameter value based on block name and LHA ID.
     *
     * This overload requires the ParameterProvider to have been constructed
     * with a specific ParameterType. The internal type is then combined with
     * the (block, id) to build a temporary ParamId.
     *
     * @param block  The name of the parameter block.
     * @param id     The LHA ID of the parameter.
     * @param d_type Type of data requested (central value, stat error, syst error, etc.).
     * @return The requested scalar value.
     */
    scalar_t operator()(const std::string& block, const LhaID& id, DataType d_type=DataType::VALUE) const;

    /**
     * @brief Checks if a parameter identified by ParamId exists.
     *
     * @param pid The parameter ID.
     * @return True if the parameter exists, false otherwise.
     */
    bool exists(const ParamId& pid) const;

    /**
     * @brief Checks if a parameter identified by ParamId exists.
     *
     * @param pid The parameter ID.
     * @return True if the parameter exists, false otherwise.
     */
    bool exists(const std::string& block, const LhaID& id) const;

    /**
     * @brief Checks if a parameter identified by block name and LHA ID exists.
     *
     * This call uses the internal ParameterType that was provided at
     * construction time.
     *
     * @param block The name of the parameter block.
     * @param id    The LHA ID of the parameter.
     * @return True if the parameter exists, false otherwise.
     */
    double get_scale(const std::string& block) const;

    /**
     * @brief Retrieves the renormalization scale associated with a given block.
     *
     * Delegates to Parameters::get_block_scale for the provider's ParameterType.
     *
     * @param block The name of the block.
     * @return The scale of the block.
     */
    ParameterType get_type() const;

    /**
     * @brief Retrieves the actual Parameter object corresponding to the given ParamId.
     *
     * @param pid The parameter ID.
     * @return A shared pointer to the Parameter object.
     */
    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const;

private:
    /// Optional type handled by this provider (required for block-based overloads).
    std::optional<ParameterType> p_type;

    /**
     * @brief Retrieves the value corresponding to the given ParamId and data type.
     *
     * Helper that dispatches between:
     *  - central value,
     *  - statistical uncertainty,
     *  - systematic uncertainty,
     *  - combined uncertainty.
     *
     * @param pid    The parameter ID.
     * @param d_type The data type requested.
     * @return The requested scalar value.
     */
    scalar_t get_value(const ParamId& pid, DataType d_type) const;
};


#endif // PARAMETERPROVIDER_H
