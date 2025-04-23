#ifndef __PARAMETERPROVIDER_H__
#define __PARAMETERPROVIDER_H__

#include "IDataProvider.h"
#include "IMonitor.h"
#include "General.h"
#include "Parameters.h"

class ParameterProvider : public IDataProvider<ParameterProvider> {
private:
    std::optional<ParameterType> p_type;

public:
    enum class DataType { VALUE, STD_STAT, STD_SYST, STD_COMBINED };

    ParameterProvider() = default;
    inline ParameterProvider(ParameterType p_type) : p_type(p_type) { Parameters::GetInstance(p_type); }

    scalar_t operator()(const ParamId& pid, DataType d_type=DataType::VALUE);
    scalar_t operator()(const std::string& block, const LhaID& id, DataType d_type=DataType::VALUE) const;

    bool exists(const ParamId& pid) const;
    bool exists(const std::string& block, const LhaID& id) const;

    ParameterType get_type() const;
    std::shared_ptr<Parameter> get_parameter(const ParamId& pid) const;
};


#endif // __PARAMETERPROVIDER_H__
