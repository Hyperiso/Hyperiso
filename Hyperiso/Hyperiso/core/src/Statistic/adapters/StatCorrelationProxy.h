#ifndef STAT_CORRELATION_PROXY_H
#define STAT_CORRELATION_PROXY_H

#include "IStatParameterProxy.h"
#include "ParameterProvider.h"
#include "CorrelationProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"

class StatCorrelationProxy : public IStatParameterProxy<std::string, LhaID> {
public:
    double operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationProvider::CorrelationType type);
    double operator()(const Observables& pid_1, const Observables& pid_2, CorrelationProvider::CorrelationType type);
    double operator()(const ObservableId& pid_1, const ObservableId& pid_2, CorrelationProvider::CorrelationType type);

    
private:
    CorrelationProvider cp;
};

#endif 