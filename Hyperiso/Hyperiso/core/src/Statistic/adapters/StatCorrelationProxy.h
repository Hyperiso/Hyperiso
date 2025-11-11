#ifndef STAT_CORRELATION_PROXY_H
#define STAT_CORRELATION_PROXY_H

#include "IStatCorrelationProxy.h"
#include "CorrelationProvider.h"
#include "Include.h"

class StatCorrelationProxy  : public IStatCorrelationProxy
{
public:
    using Type = CorrelationProvider::CorrelationType;

    double operator()(const ParamId&, const ParamId&, Type) override;
    double operator()(const Observables&, const Observables&, Type) override;
    double operator()(const ObservableId&, const ObservableId&, Type) override;

private:
    CorrelationProvider cp;
};


#endif 