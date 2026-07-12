#ifndef MC_ENGINE_H
#define MC_ENGINE_H

#include <vector>
#include <random>
#include <string>
#include <memory>

#include "INuisanceSampler.h"
#include "ports/IModel.h"
#include "Statistics.h"
#include "GaussianApprox.h"
#include "Indexing.h"
#include "StatisticProgress.h"

/**
 * @file MCEngine.h
 * @brief Monte Carlo propagation of nuisance-parameter uncertainties.
 *
 * The Monte Carlo engine samples nuisance parameters, evaluates model
 * predictions, summarizes the resulting observable distributions, and builds a
 * regularized observable covariance matrix.
 *
 * @see INuisanceSampler
 * @see IModel
 * @see GaussianSummary
 */

/**
 * @struct MCConfig
 * @brief Runtime configuration for Monte Carlo nuisance propagation.
 */
struct MCConfig {
    /** Number of accepted Monte Carlo predictions to generate. */
    std::size_t draws = 10000;

    /** Absolute skewness threshold below which a distribution is treated as symmetric. */
    double skew_abs_threshold = 0.2;

    /** Relative diagonal ridge added before covariance inversion. */
    double covariance_ridge_rel = 1e-8;

    /** Absolute diagonal ridge added before covariance inversion. */
    double covariance_ridge_abs = 1e-12;

    /** Whether failed model evaluations should be rejected and retried. */
    bool retry_failed_predictions = true;

    /** Maximum number of rejected model evaluations before aborting. */
    std::size_t max_prediction_failures = 20000;

    /** Number of MC worker threads. A value <= 1 keeps the serial path. */
    std::size_t n_threads = 1;

    /** Force thread-configurable decays to @ref forced_decay_threads during parallel MC. */
    bool force_decay_threads_to_one = true;

    /** Decay thread count used while parallel MC workers are active. */
    std::size_t forced_decay_threads = 1;

    /** Print a progress line with measured draw rate and ETA. */
    bool print_progress = false;

    /** Number of first accepted draws used before the first stable ETA estimate. */
    std::size_t progress_probe_draws = 5;

    /** Accepted-draw stride between progress updates. */
    std::size_t progress_update_every = 1;

    /** Optional frontend monitor receiving thread-safe progress snapshots. */
    std::shared_ptr<StatisticProgressMonitor> progress_monitor {};

    /** Whether accepted observable samples should be written to CSV. */
    bool write_samples_csv = false;

    /** Output path used when @ref write_samples_csv is true. */
    std::string samples_csv_path = "obs_samples.csv";
};

/**
 * @struct MCRealization
 * @brief Raw Monte Carlo samples accepted by the engine.
 */
struct MCRealization {
    /** Sampled observable values, one map per accepted draw. */
    ObsSamples sampled_obss;

    /** Nuisance-parameter values used for each accepted draw. */
    NuisanceSamples sampled_params;
};

/**
 * @struct MCObservableCovariance
 * @brief Empirical observable covariance and its inverse.
 */
struct MCObservableCovariance {
    /** Observable identifiers defining the ordering of matrix rows and columns. */
    std::vector<BinnedObservableId> ids;

    /** Empirical mean for each observable in @ref ids order. */
    std::vector<double> mean;

    /** Regularized empirical covariance matrix. */
    RealMatrix covariance;

    /** Inverse of @ref covariance. */
    RealMatrix covariance_inv;
};

/**
 * @brief Builds a regularized empirical covariance matrix from observable samples.
 *
 * The samples are read in the order specified by @p ids.  A diagonal ridge
 * equal to @c max(ridge_abs, ridge_rel * max(mean_variance, 1)) is added before
 * inversion to improve numerical stability.
 *
 * @param S Monte Carlo observable samples.
 * @param ids Observable identifiers defining the covariance ordering.
 * @param ridge_rel Relative ridge factor.
 * @param ridge_abs Absolute ridge floor.
 * @return Empirical covariance, inverse covariance and column means.
 *
 * @throws std::invalid_argument if the sample set is empty, if @p ids is empty,
 *         or if fewer than two samples are provided.
 */
MCObservableCovariance covariance_from_obs_samples(
    const ObsSamples& S,
    const std::vector<BinnedObservableId>& ids,
    double ridge_rel = 1e-8,
    double ridge_abs = 1e-12
);

/**
 * @brief Extracts the observable ordering from the first Monte Carlo sample.
 *
 * @param S Monte Carlo observable samples.
 * @return Observable identifiers present in the first sample.
 *
 * @throws std::invalid_argument if @p S is empty.
 */
std::vector<BinnedObservableId> covariance_ids_from_first_sample(
    const ObsSamples& S
);

/**
 * @struct MCResult
 * @brief Complete output of a Monte Carlo propagation run.
 */
struct MCResult {
    /** Raw accepted Monte Carlo predictions and nuisance samples. */
    MCRealization mc_real;

    /** Gaussian or split-Gaussian summary of each observable distribution. */
    std::vector<GaussianSummary> summary;

    /** Empirical observable covariance matrix and inverse. */
    MCObservableCovariance covariance;
};

/**
 * @class MonteCarloEngine
 * @brief Samples nuisance parameters and propagates them through a model.
 *
 * The engine repeatedly draws nuisance values from an @ref INuisanceSampler,
 * calls @ref IModel::predict_optimized, rejects failed or non-finite
 * predictions when configured to do so, and returns both raw and summarized
 * Monte Carlo outputs.
 */
class MonteCarloEngine {
public:
    /**
     * @brief Constructs the Monte Carlo engine.
     *
     * @param model Model used to compute predictions for each nuisance draw.
     * @param sampler Nuisance sampler used to generate random parameter maps.
     * @param cfg Monte Carlo configuration.
     */
    MonteCarloEngine(const std::shared_ptr<IModel>& model, const INuisanceSampler& sampler, MCConfig cfg)
    : model_(model), sampler_(sampler), cfg_(cfg) {}

    /**
     * @brief Generates accepted model predictions for a fixed fit-parameter point.
     *
     * @param p Fit-parameter values held fixed during the Monte Carlo run.
     * @return Raw accepted observable and nuisance samples.
     *
     * @throws std::exception Re-throws model or sampler failures once the retry
     *         policy is exhausted.
     */
    MCRealization sample_predictions(const std::map<ParamId, double>& p) const;

    MCRealization sample_predictions_serial(const std::map<ParamId, double>& p) const;
    MCRealization sample_predictions_parallel(const std::map<ParamId, double>& p) const;

    /**
     * @brief Runs Monte Carlo propagation and computes summary statistics.
     *
     * In addition to returning the result, the current implementation writes the
     * accepted observable samples to @c obs_samples.csv.
     *
     * @param p Fit-parameter values held fixed during the Monte Carlo run.
     * @return Raw samples, Gaussian summaries and empirical covariance.
     */
    MCResult summarize(const std::map<ParamId, double>& p) const;

private:
    const std::shared_ptr<IModel>& model_;  ///< Model evaluated for each accepted draw.
    const INuisanceSampler& sampler_;       ///< Source of nuisance-parameter samples.
    MCConfig cfg_;                          ///< Monte Carlo configuration.
};

#endif