#ifndef __CORRELATIONPROVIDER_H__
#define __CORRELATIONPROVIDER_H__

#include "IDataProvider.h"
#include "General.h"
#include "Parameter.h"

/**
 * @class CorrelationProvider
 * @ingroup DataProvidersModule
 * @brief Provides access to statistical, systematic, or combined correlations between parameters or observables.
 */
class CorrelationProvider : public IDataProvider<CorrelationProvider> {
public:
    enum class CorrelationType { STAT, SYST, COMBINED };

    /**
     * @brief Retrieves a correlation between two parameters.
     * @param pid_1 First parameter ID.
     * @param pid_2 Second parameter ID.
     * @param type Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationType type);

    /**
     * @brief Retrieves a correlation between two observables.
     * @param pid_1 First observable.
     * @param pid_2 Second observable.
     * @param type Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const Observables& pid_1, const Observables& pid_2, CorrelationType type);
};

#endif // __CORRELATIONPROVIDER_H__
