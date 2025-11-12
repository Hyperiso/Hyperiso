#ifndef ISTAT_CORRELATION_PROXY_H
#define ISTAT_CORRELATION_PROXY_H


#include "CorrelationProvider.h"
#include "Include.h"

struct IStatCorrelationProxy {
    using Type = CorrelationProvider::CorrelationType;
    virtual ~IStatCorrelationProxy() = default;

    virtual double operator()(const ParamId&, const ParamId&, Type) = 0;
    virtual double operator()(const Observables&, const Observables&, Type) = 0;
    virtual double operator()(const ObservableId&, const ObservableId&, Type) = 0;
};

#endif