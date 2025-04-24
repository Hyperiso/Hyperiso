#ifndef CORRELATION_PROXY_H
#define CORRELATION_PROXY_H

#include "IObsParameterProxy.h"
#include "ParameterProvider.h"
#include "CorrelationProvider.h"
#include "HyperisoMaster.h"
#include "Include.h"

class CorrelationProxy : public IObsParameterProxy<std::string, LhaID> {
public:
    double operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationProvider::CorrelationType type);
    double operator()(const Observables& pid_1, const Observables& pid_2, CorrelationProvider::CorrelationType type);
    
private:
    CorrelationProvider cp;
};

#endif 