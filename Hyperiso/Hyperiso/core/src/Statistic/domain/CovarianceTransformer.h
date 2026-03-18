#ifndef COVARIANCE_TRANSFORMER_H
#define COVARIANCE_TRANSFORMER_H

#include "IStatParameterProxy.h"
#include "ParamID.h"
#include "CorrelationProvider.h"
#include "IStatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "ExperimentObs.h"

/**
 * @file CovarianceTransformer.h
 * @brief Helper used to transform collections of identifiers into correlation matrices/maps.
 *
 * This class provides convenience utilities for building dense or associative
 * correlation representations from collections of:
 * - parameters (@ref ParamId),
 * - experiment-scoped observables (@ref ExperimentObs),
 * - observables identified by @ref BinnedObservableId together with an experiment name.
 *
 * Internally, all correlation coefficients are queried through an
 * @ref IStatCorrelationProxy using
 * @ref CorrelationProvider::CorrelationType::COMBINED.
 *
 * For parameter collections, the output is built directly from the parameter IDs.
 *
 * For observable collections, two usage modes are supported:
 * - explicit mode using @ref ExperimentObs,
 * - convenience mode for a homogeneous collection of @ref BinnedObservableId
 *   belonging to the same experiment, by passing the experiment name separately.
 *
 * The class also provides a filtering helper `check_if_corr(...)` for parameters,
 * which keeps only parameters with a strictly positive combined standard
 * uncertainty as returned by the parameter proxy.
 */

/**
 * @class CovarianceTransformer
 * @brief Builds correlation tables/matrices from parameter or observable collections.
 *
 * This transformer is a thin utility layer used by the statistics code to
 * assemble correlation structures suitable for downstream covariance/copula
 * construction.
 *
 * It does not store any correlation data itself; instead, it delegates all
 * lookups to:
 * - an @ref IStatCorrelationProxy for correlations,
 * - an @ref IStatParameterProxy for parameter uncertainty checks.
 */
class CovarianceTransformer {
public:
    /**
     * @brief Constructs a covariance transformer.
     *
     * @param corr_proxy Proxy used to retrieve correlation coefficients.
     * @param par_proxy  Proxy used to retrieve parameter uncertainties.
     */
    CovarianceTransformer(std::shared_ptr<IStatCorrelationProxy> corr_proxy,
                          std::shared_ptr<IStatParameterProxy> par_proxy)
        : corr_proxy(std::move(corr_proxy)), par_proxy(std::move(par_proxy)) { }

    /**
     * @brief Builds a dense correlation matrix for a list of parameters.
     *
     * The output matrix is ordered according to the input vector. Entry `(i,j)`
     * contains the combined correlation between `ids[i]` and `ids[j]`.
     *
     * @param ids Ordered list of parameter identifiers.
     * @return Dense square matrix of combined correlations.
     */
    std::vector<std::vector<double>> transform(const std::vector<ParamId>& ids);

    /**
     * @brief Builds an associative correlation map for a set of parameters.
     *
     * The values of the input map are ignored; only the keys are used to define
     * the parameter set over which the correlation structure is constructed.
     *
     * @param ids Map whose keys define the parameter set.
     * @return Nested map of combined correlations.
     */
    std::map<ParamId, std::map<ParamId, double>> transform(const std::map<ParamId, double>& ids);

    /**
     * @brief Builds a dense correlation matrix for a list of fully explicit
     *        experiment-scoped observables.
     *
     * The output matrix is ordered according to the input vector. Entry `(i,j)`
     * contains the combined correlation between `ids[i]` and `ids[j]`.
     *
     * @param ids Ordered list of experiment-scoped observable identifiers.
     * @return Dense square matrix of combined correlations.
     */
    std::vector<std::vector<double>> transform(const std::vector<ExperimentObs>& ids);

    /**
     * @brief Builds an associative correlation map for a set of fully explicit
     *        experiment-scoped observables.
     *
     * The values of the input map are ignored; only the keys are used to define
     * the observable set over which the correlation structure is constructed.
     *
     * @param ids Map whose keys define the experiment-scoped observable set.
     * @return Nested map of combined correlations.
     */
    std::map<ExperimentObs, std::map<ExperimentObs, double>> transform(const std::map<ExperimentObs, double>& ids);

    /**
     * @brief Builds an associative correlation map for a set of binned observables
     *        belonging to the same experiment.
     *
     * This is a convenience overload around the @ref ExperimentObs-based API.
     *
     * @param experiment Experiment name used to scope all observables.
     * @param ids Map whose keys define the observable set.
     * @return Nested map of combined correlations.
     */
    std::map<BinnedObservableId, std::map<BinnedObservableId, double>>
    transform(const std::string& experiment, const std::map<BinnedObservableId, double>& ids);

    /**
     * @brief Builds a dense correlation matrix for a list of binned observables
     *        belonging to the same experiment.
     *
     * This is a convenience overload around the @ref ExperimentObs-based API.
     *
     * @param experiment Experiment name used to scope all observables.
     * @param ids Ordered list of binned observable identifiers.
     * @return Dense square matrix of combined correlations.
     */
    std::vector<std::vector<double>>
    transform(const std::string& experiment, const std::vector<BinnedObservableId>& ids);

    /**
     * @brief Filters a parameter list, keeping only entries with positive
     *        combined uncertainty.
     *
     * A parameter is kept if `(*par_proxy)(elem, DataType::STD_COMBINED) > 0`.
     *
     * @param ids Input parameter list.
     * @return Filtered list of parameters considered correlated / relevant.
     */
    std::vector<ParamId> check_if_corr(const std::vector<ParamId>& ids);

private:
    std::shared_ptr<IStatCorrelationProxy> corr_proxy; ///< Proxy used to query correlations.
    std::shared_ptr<IStatParameterProxy> par_proxy;    ///< Proxy used to query parameter uncertainties.
};

#endif // COVARIANCE_TRANSFORMER_H