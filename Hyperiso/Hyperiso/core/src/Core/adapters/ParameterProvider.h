#ifndef __PARAMETERPROVIDER_H__
#define __PARAMETERPROVIDER_H__

#include "IDataProvider.h"
#include "IMonitor.h"
#include "General.h"

class ParameterProvider : public IDataProvider<ParameterProvider> {
private:
    ParameterType* p_type {nullptr};

public:
    enum class DataType { VALUE, STD_STAT, STD_SYST, STD_COMBINED };

    ParameterProvider() = default;
    inline ParameterProvider(ParameterType p_type) : p_type(&p_type) {}

    double operator()(const ParamId& pid, DataType d_type=DataType::VALUE);
    double operator()(const std::string& block, const LhaID& id, DataType d_type=DataType::VALUE);
};


#endif // __PARAMETERPROVIDER_H__
