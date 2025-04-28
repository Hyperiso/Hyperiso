#ifndef __PARAMETERPROVIDER_H__
#define __PARAMETERPROVIDER_H__

#include "IDataProvider.h"
#include "IMonitor.h"
#include "General.h"
#include "Parameters.h"

/**
 * @class ParameterProvider
 * @ingroup DataProvidersModule
 * @brief Provides access to parameter values, errors, and existence checks.
 */
class ParameterProvider : public IDataProvider<ParameterProvider> {
public:
    enum class DataType { VALUE, STD_STAT, STD_SYST, STD_COMBINED };

    ParameterProvider() = default;
    inline ParameterProvider(ParameterType p_type) : p_type(p_type) { Parameters::GetInstance(p_type); }

    /**
     * @brief Retrieves a parameter value based on ParamId and requested data type.
     * @param pid The parameter ID.
     * @param d_type Type of data requested (value, stat error, syst error, etc.).
     * @return The requested scalar value.
     */
    scalar_t operator()(const ParamId& pid, DataType d_type=DataType::VALUE);

    /**
     * @brief Retrieves a parameter value based on block name and LHA ID.
     * @param block The name of the parameter block.
     * @param id The LHA ID of the parameter.
     * @param d_type Type of data requested (value, stat error, syst error, etc.).
     * @return The requested scalar value.
     */
    scalar_t operator()(const std::string& block, const LhaID& id, DataType d_type=DataType::VALUE) const;

    /**
     * @brief Checks if a parameter identified by ParamId exists.
     * @param pid The parameter ID.
     * @return True if the parameter exists, false otherwise.
     */
    bool exists(const ParamId& pid) const;

    /**
     * @brief Checks if a parameter identified by block name and LHA ID exists.
     * @param block The name of the parameter block.
     * @param id The LHA ID of the parameter.
     * @return True if the parameter exists, false otherwise.
     */
    bool exists(const std::string& block, const LhaID& id) const;

    /**
     * @brief Retrieves the ParameterType managed by this provider.
     * @return The ParameterType (e.g., SM, BSM, etc.).
     */
    ParameterType get_type() const;

    /**
     * @brief Retrieves the actual Parameter object corresponding to the given ParamId.
     * @param pid The parameter ID.
     * @return A shared pointer to the Parameter object.
     */
    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const;

private:
    std::optional<ParameterType> p_type;

    /**
     * @brief Retrieves the value corresponding to the given ParamId and data type (Value, std_stat/syst/combined).
     * @param pid The parameter ID.
     * @param d_type The data type requested.
     * @return The requested scalar value.
     */
    scalar_t get_value(const ParamId& pid, DataType d_type) const;
};


#endif // __PARAMETERPROVIDER_H__
