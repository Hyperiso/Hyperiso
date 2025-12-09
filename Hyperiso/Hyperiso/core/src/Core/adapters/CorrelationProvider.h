#ifndef __CORRELATIONPROVIDER_H__
#define __CORRELATIONPROVIDER_H__

#include "IDataProvider.h"
#include "ParamID.h"
#include "Parameter.h"
#include "CorrelationRepo.h"
#include "MemoryManager.h"

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

    double operator()(const ObservableId& pid_1, const ObservableId& pid_2, CorrelationType type);

    /**
     * @brief Check if correlations exist between two parameters.
     * @param pid_1 First parameter ID.
     * @param pid_2 Second parameter ID.
     * @param type Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero (exist)
     */
    bool exists(const ParamId& pid_1, const ParamId& pid_2, CorrelationType type);

    /**
     * @brief Check if correlations exist between two observables.
     * @param pid_1 First observable.
     * @param pid_2 Second observable.
     * @param type Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero (exist)
     */
    bool exists(const Observables& pid_1, const Observables& pid_2, CorrelationType type);

    bool exists(const ObservableId& pid_1, const ObservableId& pid_2, CorrelationType type);

private:
    /**
     * @brief Internal utility function to compute the correlation between two entities.
     * 
     * This function is used internally to calculate the correlation value between
     * two entities of the same type, such as parameters or observables.
     * 
     * @tparam T The type of the entities (e.g., ParamId or Observables).
     * @param id_1 The first entity.
     * @param id_2 The second entity.
     * @param type The type of correlation to compute (stat, syst, combined).
     * @return The computed correlation value.
     */
    template<typename T>
    double get_correlation(const T& id_1, const T& id_2, CorrelationType type) const;
};

#endif // __CORRELATIONPROVIDER_H__
