#ifndef MODEL_PARAM_ADAPTER_H
#define MODEL_PARAM_ADAPTER_H

#include "IDataBaseAdapter.h"
#include "ParameterProvider.h"
#include "CorrelationProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"

class CorrelationProxy : public IDataBaseProxy<std::string, LhaID> {
public:
    CorrelationProxy(ParameterType type);

    double operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationProvider::CorrelationType type);
    double operator()(const Observables& pid_1, const Observables& pid_2, CorrelationProvider::CorrelationType type);
    
private:
    CorrelationProvider cp;
};

#endif 