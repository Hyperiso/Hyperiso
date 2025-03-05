#ifndef __CORRELATIONPROVIDER_H__
#define __CORRELATIONPROVIDER_H__

#include "IDataProvider.h"
#include "General.h"
#include "Parameter.h"

class CorrelationProvider : public IDataProvider<CorrelationProvider> {
public:
    enum class CorrelationType { STAT, SYST, COMBINED };

    double operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationType type);
    double operator()(const Observables& pid_1, const Observables& pid_2, CorrelationType type);
};

#endif // __CORRELATIONPROVIDER_H__
