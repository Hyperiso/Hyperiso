#ifndef CORRELATIONPROVIDER_H
#define CORRELATIONPROVIDER_H

#include "IDataProvider.h"
#include "ParamID.h"
#include "Parameter.h"
#include "CorrelationRepo.h"
#include "MemoryManager.h"

/**
 * @file CorrelationProvider.h
 * @brief High-level access to parameter and observable correlations.
 *
 * This header declares the CorrelationProvider class, a concrete
 * IDataProvider that exposes a uniform interface to:
 *  - statistical correlations,
 *  - systematic correlations,
 *  - combined correlations (hypotenuse of stat and syst),
 * between:
 *  - parameters (ParamId),
 *  - observables (Observables or ObservableId).
 *
 * Internally it relies on the CorrelationRepository owned by the
 * MemoryManager singleton.
 */

/**
 * @class CorrelationProvider
 * @ingroup DataProvidersModule
 * @brief Provides access to statistical, systematic, or combined correlations between parameters or observables.
 *
 * CorrelationProvider implements a callable interface:
 *  - `operator()(ParamId, ParamId, CorrelationType)` for parameters,
 *  - `operator()(Observables, Observables, CorrelationType)` for enum-based observables,
 *  - `operator()(ObservableId, ObservableId, CorrelationType)` for explicit observable IDs.
 *
 * It delegates the actual lookup and combination to the global
 * CorrelationRepository exposed by MemoryManager.
 *
 * Example:
 * @code
 * CorrelationProvider corr;
 * ParamId p1{ParameterType::SM, "MASS", LhaID(5)};
 * ParamId p2{ParameterType::SM, "MASS", LhaID(6)};
 *
 * double rho_stat = corr(p1, p2, CorrelationProvider::CorrelationType::STAT);
 * bool   exists   = corr.exists(p1, p2, CorrelationProvider::CorrelationType::COMBINED);
 * @endcode
 */
class CorrelationProvider : public IDataProvider<CorrelationProvider> {
public:
    /// Type of correlation component requested.
    enum class CorrelationType { STAT, SYST, COMBINED };

    /**
     * @brief Retrieves a correlation between two parameters.
     *
     * @param pid_1 First parameter ID.
     * @param pid_2 Second parameter ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const ParamId& pid_1, const ParamId& pid_2, CorrelationType type);

    /**
     * @brief Retrieves a correlation between two observables specified by enum.
     *
     * This overload accepts the high-level Observables enum and internally
     * maps it to ObservableId via ObservableMapper.
     *
     * @param pid_1 First observable (enum).
     * @param pid_2 Second observable (enum).
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const Observables& pid_1, const Observables& pid_2, CorrelationType type);

    /**
     * @brief Retrieves a correlation between two observables specified by ObservableId.
     *
     * @param pid_1 First observable ID.
     * @param pid_2 Second observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const ObservableId& pid_1, const ObservableId& pid_2, CorrelationType type);

    /**
     * @brief Retrieves a correlation between two observables specified by BinnedObservableId.
     *
     * @param pid_1 First observable ID.
     * @param pid_2 Second observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const BinnedObservableId& pid_1, const BinnedObservableId& pid_2, CorrelationType type);

    /**
     * @brief Checks if a correlation between two parameters is non-zero.
     *
     * Internally calls the corresponding operator() and tests whether the
     * returned value is non-zero.
     *
     * @param pid_1 First parameter ID.
     * @param pid_2 Second parameter ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero, false otherwise.
     */
    bool exists(const ParamId& pid_1, const ParamId& pid_2, CorrelationType type);

    /**
     * @brief Checks if a correlation between two observables (enum) is non-zero.
     *
     * @param pid_1 First observable (enum).
     * @param pid_2 Second observable (enum).
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero, false otherwise.
     */
    bool exists(const Observables& pid_1, const Observables& pid_2, CorrelationType type);

    /**
     * @brief Checks if a correlation between two observables (ObservableId) is non-zero.
     *
     * @param pid_1 First observable ID.
     * @param pid_2 Second observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero, false otherwise.
     */
    bool exists(const ObservableId& pid_1, const ObservableId& pid_2, CorrelationType type);

    /**
     * @brief Checks if a correlation between two observables (BinnedObservableId) is non-zero.
     *
     * @param pid_1 First observable ID.
     * @param pid_2 Second observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero, false otherwise.
     */
    bool exists(const BinnedObservableId& pid_1, const BinnedObservableId& pid_2, CorrelationType type);

private:
    /**
     * @brief Internal utility function to compute the correlation between two entities.
     *
     * This function is used internally to calculate the correlation value between
     * two entities of the same type, such as parameters or observables. It handles
     * the splitting between statistical, systematic, and combined correlations.
     *
     * @tparam T  The type of the entities (e.g., ParamId, Observables, ObservableId).
     * @param id_1 The first entity.
     * @param id_2 The second entity.
     * @param type The type of correlation to compute (stat, syst, combined).
     * @return The computed correlation value.
     */
    template<typename T>
    double get_correlation(const T& id_1, const T& id_2, CorrelationType type) const;
};

#endif // CORRELATIONPROVIDER_H
