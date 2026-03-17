#ifndef CORRELATIONPROVIDER_H
#define CORRELATIONPROVIDER_H

#include "IDataProvider.h"
#include "ParamID.h"
#include "Parameter.h"
#include "CorrelationRepo.h"
#include "MemoryManager.h"
#include "ExperimentObs.h"

/**
 * @file CorrelationProvider.h
 * @brief High-level access to parameter and experiment-scoped observable correlations.
 *
 * This header declares the CorrelationProvider class, a concrete
 * IDataProvider that exposes a uniform interface to:
 *  - statistical correlations,
 *  - systematic correlations,
 *  - combined correlations (hypotenuse of stat and syst),
 * between:
 *  - parameters (ParamId),
 *  - experiment-scoped observables (ExperimentObs),
 *  - observables identified by enum or ID, together with an experiment name.
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
 *  - `operator()(ExperimentObs, ExperimentObs, CorrelationType)` for fully explicit experiment-scoped observables,
 *  - `operator()(experiment, Observables, Observables, CorrelationType)` for enum-based observables,
 *  - `operator()(experiment, ObservableId, ObservableId, CorrelationType)` for explicit observable IDs,
 *  - `operator()(experiment, BinnedObservableId, BinnedObservableId, CorrelationType)` for binned observables,
 *  - `operator()(exp1, BinnedObservableId, exp2, BinnedObservableId, CorrelationType)` for cross-experiment queries.
 *
 * It delegates the actual lookup and combination to the global
 * CorrelationRepository exposed by MemoryManager.
 *
 * Example:
 * @code
 * CorrelationProvider corr;
 *
 * ParamId p1{ParameterType::SM, "MASS", LhaID(5)};
 * ParamId p2{ParameterType::SM, "MASS", LhaID(6)};
 * double rho_stat_param = corr(p1, p2, CorrelationProvider::CorrelationType::STAT);
 *
 * BinnedObservableId o1 = ...;
 * BinnedObservableId o2 = ...;
 * double rho_stat_obs = corr("LHCb", o1, o2, CorrelationProvider::CorrelationType::STAT);
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
     * @brief Retrieves a correlation between two fully explicit experiment-scoped observables.
     *
     * @param oid_1 First experiment-scoped observable.
     * @param oid_2 Second experiment-scoped observable.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const ExperimentObs& oid_1, const ExperimentObs& oid_2, CorrelationType type);

    /**
     * @brief Retrieves a correlation between two observables specified by enum,
     *        within the same experiment.
     *
     * This overload accepts the high-level Observables enum and internally
     * maps it to ObservableId via ObservableMapper.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param oid_1 First observable (enum).
     * @param oid_2 Second observable (enum).
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const std::string& experiment,
                      const Observables& oid_1,
                      const Observables& oid_2,
                      CorrelationType type);

    /**
     * @brief Retrieves a correlation between two observables specified by ObservableId,
     *        within the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param oid_1 First observable ID.
     * @param oid_2 Second observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const std::string& experiment,
                      const ObservableId& oid_1,
                      const ObservableId& oid_2,
                      CorrelationType type);

    /**
     * @brief Retrieves a correlation between two observables specified by BinnedObservableId,
     *        within the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param oid_1 First binned observable ID.
     * @param oid_2 Second binned observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const std::string& experiment,
                      const BinnedObservableId& oid_1,
                      const BinnedObservableId& oid_2,
                      CorrelationType type);

    /**
     * @brief Retrieves a correlation between two binned observables coming
     *        from possibly different experiments.
     *
     * @param exp_1 Experiment associated with the first observable.
     * @param oid_1 First binned observable ID.
     * @param exp_2 Experiment associated with the second observable.
     * @param oid_2 Second binned observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return The requested correlation value.
     */
    double operator()(const std::string& exp_1,
                      const BinnedObservableId& oid_1,
                      const std::string& exp_2,
                      const BinnedObservableId& oid_2,
                      CorrelationType type);

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
     * @brief Checks if a correlation between two experiment-scoped observables is non-zero.
     *
     * @param oid_1 First experiment-scoped observable.
     * @param oid_2 Second experiment-scoped observable.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero, false otherwise.
     */
    bool exists(const ExperimentObs& oid_1, const ExperimentObs& oid_2, CorrelationType type);

    /**
     * @brief Checks if a correlation between two observables (enum) in the same experiment is non-zero.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param oid_1 First observable (enum).
     * @param oid_2 Second observable (enum).
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero, false otherwise.
     */
    bool exists(const std::string& experiment,
                const Observables& oid_1,
                const Observables& oid_2,
                CorrelationType type);

    /**
     * @brief Checks if a correlation between two observables (ObservableId) in the same experiment is non-zero.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param oid_1 First observable ID.
     * @param oid_2 Second observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero, false otherwise.
     */
    bool exists(const std::string& experiment,
                const ObservableId& oid_1,
                const ObservableId& oid_2,
                CorrelationType type);

    /**
     * @brief Checks if a correlation between two observables (BinnedObservableId)
     *        in the same experiment is non-zero.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param oid_1 First binned observable ID.
     * @param oid_2 Second binned observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero, false otherwise.
     */
    bool exists(const std::string& experiment,
                const BinnedObservableId& oid_1,
                const BinnedObservableId& oid_2,
                CorrelationType type);

    /**
     * @brief Checks if a correlation between two binned observables coming from
     *        possibly different experiments is non-zero.
     *
     * @param exp_1 Experiment associated with the first observable.
     * @param oid_1 First binned observable ID.
     * @param exp_2 Experiment associated with the second observable.
     * @param oid_2 Second binned observable ID.
     * @param type  Type of correlation requested (stat, syst, combined).
     * @return True if the correlation is non-zero, false otherwise.
     */
    bool exists(const std::string& exp_1,
                const BinnedObservableId& oid_1,
                const std::string& exp_2,
                const BinnedObservableId& oid_2,
                CorrelationType type);

private:
    /**
     * @brief Internal utility function to compute the correlation between two entities.
     *
     * This helper is used for homogeneous key types directly understood by
     * CorrelationRepository, namely:
     *   - ParamId
     *   - ExperimentObs
     *
     * It handles the splitting between statistical, systematic, and combined
     * correlations.
     *
     * @tparam T  The type of the entities.
     * @param id_1 The first entity.
     * @param id_2 The second entity.
     * @param type The type of correlation to compute (stat, syst, combined).
     * @return The computed correlation value.
     */
    template<typename T>
    double get_correlation(const T& id_1, const T& id_2, CorrelationType type) const;

    /**
     * @brief Internal utility function to compute the correlation between two
     *        observables in the same experiment.
     *
     * @param experiment Experiment name used to scope both observables.
     * @param id_1 First binned observable ID.
     * @param id_2 Second binned observable ID.
     * @param type Type of correlation to compute (stat, syst, combined).
     * @return The computed correlation value.
     */
    double get_correlation(const std::string& experiment,
                           const BinnedObservableId& id_1,
                           const BinnedObservableId& id_2,
                           CorrelationType type) const;

    /**
     * @brief Internal utility function to compute the correlation between two
     *        observables possibly associated with different experiments.
     *
     * @param exp_1 Experiment associated with the first observable.
     * @param id_1  First binned observable ID.
     * @param exp_2 Experiment associated with the second observable.
     * @param id_2  Second binned observable ID.
     * @param type  Type of correlation to compute (stat, syst, combined).
     * @return The computed correlation value.
     */
    double get_correlation(const std::string& exp_1,
                           const BinnedObservableId& id_1,
                           const std::string& exp_2,
                           const BinnedObservableId& id_2,
                           CorrelationType type) const;
};

#endif // CORRELATIONPROVIDER_H