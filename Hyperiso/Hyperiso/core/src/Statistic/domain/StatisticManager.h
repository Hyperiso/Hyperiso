#ifndef STATISTIC_MANAGER_H
#define STATISTIC_MANAGER_H

#include <vector>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <optional>
#include <limits>
#include <numeric>
#include <set>
#include <random>
#include <map>

#include "Include.h"
#include "IModel.h"
#include "IStatCorrelationProxy.h"
#include "IStatParameterProxy.h"
#include "IStatSourcesProxy.h"
#include "IStatDependencyPruner.h"
#include "CovarianceTransformer.h"
#include "JointDistribution.h"
#include "RvgNuisanceSampler.h"
#include "MCEngine.h"
#include "MarginalFactory.h"
#include "CopulaFactory.h"
#include "Fit.h"
#include "ChiSquaredLikelihood.h"
#include "MarginalConfigFactory.h"
#include "INuisanceReader.h"
#include "NuisanceSpec.h"
#include "IStatParamOptimizerProxy.h"
#include "ChiSquaredLikelihood.h"

/**
 * @file StatisticManager.h
 * @brief High-level orchestration of statistical uncertainty propagation, likelihood construction and fit scans.
 *
 * The statistic manager is the central bridge between the physics model, statistical
 * proxies, nuisance-parameter definitions and the fitting/contouring machinery. It
 * builds the cached observable and nuisance state, constructs likelihood objects,
 * runs maximum-likelihood fits, computes Monte Carlo uncertainty summaries and
 * exposes confidence-contour and likelihood-scan utilities.
 *
 * @see MLFitter
 * @see MonteCarloEngine
 * @see JointDistribution
 * @see NuisanceSpec
 */

/**
 * @enum StatisticLikelihoodMode
 * @brief Selects the likelihood backend used by StatisticManager::compute_MLE().
 */
enum class StatisticLikelihoodMode {
    PROFILED_NUISANCE,  ///< Full likelihood with explicit nuisance parameters profiled during the fit.
    CHI2_MC_COVARIANCE  ///< Fast chi-square likelihood using MC theory covariance plus experimental covariance.
};

/**
 * @struct StatisticConfig
 * @brief Runtime configuration for statistical propagation, fitting and nuisance handling.
 *
 * The defaults are intentionally conservative: a Gaussian copula is used for both
 * nuisance and experimental distributions, HESSE is enabled for the MLE, and local
 * sensitivity pruning is enabled to reduce the nuisance set before expensive fits.
 */
struct AdvancedStatisticConfig {
    std::map<ParamId, MarginalType> override_nuisance_marginals {};      ///< Per-parameter overrides for nuisance marginal laws.
    std::map<ExperimentObs, MarginalType> override_exp_data_marginals {};///< Per-observable overrides for experimental-data marginals.
    CopulaType nuisance_copula_type = CopulaType::GAUSSIAN;              ///< Copula used to correlate nuisance parameters.
    CopulaType exp_data_copula_type = CopulaType::GAUSSIAN;              ///< Copula used to correlate experimental observables.

    std::size_t MLE_max_iter = 500;                                      ///< Maximum number of minimizer function calls/iterations.
    double MLE_tol = 1e-8;                                               ///< Minimizer tolerance passed to the backend.
    unsigned MLE_strategy = 2;                                           ///< Backend minimization strategy; zero means backend default where supported.
    bool MLE_run_hesse = true;                                           ///< Whether to request HESSE/covariance estimation after the fit.
    bool MLE_request_minos = false;                                      ///< Whether to request MINOS errors when supported by the backend build.
    bool MLE_verbose = false;                                            ///< Enables verbose output from the fit backend.

    double nuisance_relevance_cutoff = 1e-8;                             ///< Relative-uncertainty cutoff for the first nuisance preselection pass.
    bool nuisance_sensitivity_pruning = true;                            ///< Enables local model-sensitivity pruning of nuisance candidates.
    double nuisance_sensitivity_probe_sigmas = 1.0;                      ///< Size of the +/- finite-difference probe in units of nuisance sigma.
    double nuisance_sensitivity_rel_cutoff = 1e-6;                       ///< Relative observable shift required to keep a nuisance.
    double nuisance_sensitivity_abs_cutoff = 1e-12;                      ///< Absolute observable shift required to keep a nuisance.
    double nuisance_sensitivity_scale_floor = 1e-3;                      ///< Lower scale used when normalizing relative observable shifts.
    int nuisance_sensitivity_contexts = 2;                               ///< Number of contexts tested by sensitivity pruning; negative disables the check.
    double nuisance_sensitivity_context_sigma = 0.35;                    ///< Randomized-context spread in nuisance sigma units.
    unsigned nuisance_sensitivity_seed = 12345;                          ///< RNG seed used to build sensitivity-pruning contexts.
    bool nuisance_sensitivity_keep_on_failure = true;                    ///< Keeps a nuisance if its sensitivity probe fails.

    bool MLE_trace_first_evals = false;                                  ///< Enables debug tracing of the first likelihood evaluations.
    std::size_t MLE_trace_max_evals = 25;                                ///< Maximum number of likelihood evaluations printed when tracing is enabled.
    bool MLE_allow_profile_hessian_fallback = true;                      ///< Allows numerical profile-Hessian covariance fallback if backend covariance fails.
    double MLE_profile_hessian_step_scale = 1.0;                         ///< Step scaling used by the numerical profile-Hessian fallback.
    double MLE_profile_hessian_eig_floor_rel = 1e-8;                     ///< Relative eigenvalue floor used to regularize the fallback Hessian.

    StatisticLikelihoodMode likelihood_mode = StatisticLikelihoodMode::PROFILED_NUISANCE; ///< Likelihood mode used by compute_MLE().
    double chi2_covariance_ridge_rel = 1e-8;                             ///< Relative diagonal ridge used before inverting chi-square covariance matrices.
    double chi2_covariance_ridge_abs = 1e-12;                            ///< Absolute diagonal ridge used before inverting chi-square covariance matrices.

    bool MC_force_decay_threads_to_one = true;                           ///< Give MC priority over internal decay parallelism.
    std::size_t MC_forced_decay_threads = 1;                             ///< Decay thread count while MC workers are running.
};

/**
 * @struct StatisticConfig
 * @brief Basic runtime configuration for statistical propagation.
 *
 * Keep this structure small and user-facing: it contains the knobs that are
 * commonly changed in scripts and examples.  Advanced fit/pruning/covariance
 * controls live in @ref AdvancedStatisticConfig and are grouped under
 * @ref advanced to avoid an overloaded top-level config object.
 */
struct StatisticConfig {
    std::size_t MC_draws = 100;                 ///< Number of accepted MC draws used for uncertainty propagation.
    std::size_t MC_threads = 1;                 ///< Number of worker threads used by MC propagation.
    unsigned int MC_seed = 123456u;             ///< RNG seed used for reproducible MC nuisance and experimental-data sampling.
    double skew_abs_threshold = 0.2;            ///< Absolute skewness threshold below which a summary is treated as symmetric.

    bool print_mc_progress = true;              ///< Print MC progress with ETA based on measured draw time.
    bool print_chi2_pipeline_progress = false;  ///< Print chi-square workflow stages after/beside the MC progress bar.
    bool print_mc_config = false;               ///< Print nuisance candidates and retained MC marginal configuration.
    bool print_fit_summary = false;             ///< Print high-level fit backend summaries.
    bool print_scan_summary = false;            ///< Print likelihood-scan summaries.
    bool print_cache_summary = false;           ///< Print internal cache diagnostics.
    bool print_debug = false;                   ///< Master debug flag for low-level diagnostic output.

    bool write_mc_samples_csv = false;          ///< Write accepted MC observable samples to CSV.
    std::string mc_samples_csv_path = "obs_samples.csv"; ///< Output CSV path used when @ref write_mc_samples_csv is true.
    std::size_t mc_progress_probe_draws = 5;    ///< Number of first accepted draws used to stabilize the first ETA.
    std::size_t mc_progress_update_every = 1;   ///< Accepted-draw stride between progress updates.

    std::shared_ptr<StatisticProgressMonitor> progress_monitor {}; ///< Optional thread-safe progress sink for GUI/notebook frontends.

    std::map<ParamId, std::pair<double, double>> fit_parameter_bounds {}; ///< Optional explicit minimizer bounds keyed by fit ParamId.
    std::map<ParamId, double> fit_parameter_offsets {}; ///< Optional affine display offsets: model value = fitted value - offset.

    AdvancedStatisticConfig advanced {};        ///< Advanced fit/pruning/covariance configuration.
};

/**
 * @struct LikelihoodScanPoint
 * @brief Single point of a two-dimensional likelihood scan.
 */
struct LikelihoodScanPoint {
    double x = 0.0;             ///< Coordinate along the first scanned fit parameter.
    double y = 0.0;             ///< Coordinate along the second scanned fit parameter.
    double nll = 0.0;           ///< Negative log-likelihood value at this point.
    double delta_nll = 0.0;     ///< Difference between @ref nll and the minimum NLL found on the scan grid.
};

/**
 * @struct LikelihoodScanGrid
 * @brief Regular two-dimensional grid of likelihood-scan evaluations.
 */
struct LikelihoodScanGrid {
    ParamId x_param;                            ///< Identifier of the first scanned parameter.
    ParamId y_param;                            ///< Identifier of the second scanned parameter.
    double x_center = 0.0;                      ///< Reference value used as the scan center on x.
    double y_center = 0.0;                      ///< Reference value used as the scan center on y.
    std::size_t nx = 0;                         ///< Number of grid points along x.
    std::size_t ny = 0;                         ///< Number of grid points along y.
    std::vector<LikelihoodScanPoint> points;    ///< Flattened grid points in x-major order.
};

/**
 * @struct FitResultWithMaps
 * @brief User-facing MLE result keyed by physics parameter identifiers.
 */
struct FitResultWithMaps {
    bool fit_ok {false};                                             ///< True when the fit returned a usable parameter estimate.
    std::map<ParamId, double> p_hat;                                 ///< Best-fit values of fitted parameters.
    std::map<ParamId, double> eta_hat;                               ///< Best-fit/profiler values of nuisance parameters.
    std::map<ParamId, double> p_hat_std;                             ///< Profiled standard deviations of fitted parameters.
    std::map<ParamId, std::map<ParamId, double>> p_correlations;     ///< Correlation matrix of fitted parameters.
    double ell_hat {0.0};                                            ///< Minimum negative log-likelihood value.
};

/**
 * @struct StatCache
 * @brief Internal cache of currently selected observables, nuisances, correlations and fit state.
 */
struct StatCache {
    std::map<ExperimentObs, double> exp_obs;                         ///< Experimental central values used by the current fit.
    std::map<ParamId, double> eta_specs_real;                        ///< Selected nuisance parameters and their central values.
    std::map<ParamId, double> p_specs;                               ///< Selected fit parameters and their initial values.
    std::map<ParamId, std::map<ParamId, double>> SigmaEta;           ///< Correlation matrix of selected nuisances.
    std::map<ExperimentObs, std::map<ExperimentObs, double>> SigmaObs;///< Correlation matrix of selected experimental observables.

    FitResultWithMaps mle_result;                                    ///< Last MLE result expressed with map-based identifiers.
};

/**
 * @class StatisticManager
 * @brief Coordinates statistical inputs, nuisance distributions, MLE fits and contour/scan computations.
 *
 * The manager owns no low-level physics data itself. Instead, it delegates model
 * predictions to @ref IModel, parameter values and uncertainties to the statistical
 * proxies, and nuisance metadata to @ref INuisanceReader. Its role is to assemble
 * those sources into the distributions and likelihoods required by the fitting
 * workflow.
 */
class StatisticManager {
public:
    /**
     * @brief Constructs a statistic manager and initializes observable/nuisance state.
     *
     * The constructor computes the model observables once, loads default and user
     * nuisance specifications, merges them, and invalidates any previous fit state.
     *
     * @param config Runtime statistical configuration.
     * @param obs_int Model used to compute and predict observables.
     * @param pscp Correlation proxy for parameters and observables.
     * @param pspp Parameter proxy providing central values and uncertainties.
     * @param sp Source proxy used to resolve leaf dependencies.
     * @param dp Dependency pruner used to detach fit parameters from nuisance dependencies.
     * @param nuisance_reader Reader used to load default/user nuisance specifications.
     * @param spop Optional/statistical parameter optimizer proxy used by the surrounding workflow.
     *
     * @throws std::invalid_argument if @p nuisance_reader is null.
     */
    StatisticManager(StatisticConfig config,
                     std::shared_ptr<IModel> obs_int,
                     std::shared_ptr<IStatCorrelationProxy> pscp,
                     std::shared_ptr<IStatParameterProxy> pspp,
                     std::shared_ptr<IStatSourcesProxy> sp,
                     std::shared_ptr<IStatDependencyPruner> dp,
                     std::shared_ptr<INuisanceReader> nuisance_reader,
                     std::shared_ptr<IStatParamOptimizerProxy> spop);
    
    /**
     * @brief Builds marginal distributions for all currently cached nuisances.
     * @return Marginal distributions ordered consistently with @ref StatCache::eta_specs_real.
     */
    std::vector<std::unique_ptr<IMarginalDistribution>> build_nuisance_marginal_distributions();

    /**
     * @brief Builds the joint nuisance distribution from cached nuisance marginals and correlations.
     * @return Joint distribution over selected nuisance parameters.
     */
    std::unique_ptr<JointDistribution> build_nuisance_distribution();

    /**
     * @brief Builds the joint experimental-data distribution from cached observables and correlations.
     * @return Joint distribution over selected experimental observables.
     */
    std::unique_ptr<JointDistribution> build_exp_data_distribution();

    /**
     * @brief Reloads default and user nuisance specifications and invalidates fit state.
     */
    void reload_nuisance_specs();

    /**
     * @brief Selects a custom user nuisance-definition file and reloads nuisance specifications.
     * @param user_yaml_path Path to the user YAML/JSON nuisance file.
     */
    void set_nuisance_user_file(const fs::path& user_yaml_path);

    /**
     * @brief Clears the custom user nuisance file and reloads the default configured user file.
     */
    void clear_nuisance_user_file();

    /** @return Default nuisance specifications loaded from the default nuisance file. */
    const NuisanceRegistry& default_nuisance_specs() const { return default_nuisance_specs_; }

    /** @return User nuisance specifications loaded from the active user nuisance file. */
    const NuisanceRegistry& user_nuisance_specs() const { return user_nuisance_specs_; }

    /** @return Merged nuisance specifications, with user entries overriding default ones. */
    const NuisanceRegistry& merged_nuisance_specs() const { return merged_nuisance_specs_; }

    /**
     * @brief Computes Gaussian summaries for MC-propagated observable uncertainties.
     * @return Map from binned observable identifiers to Gaussian summaries.
     */
    std::map<BinnedObservableId, GaussianSummary> compute_uncertainties();

    /**
     * @brief Runs MC uncertainty propagation and returns both samples and summaries.
     * @return Full MC result including accepted samples, Gaussian summaries and covariance.
     */
    MCResult compute_uncertainties_and_sampling();
    
    /**
     * @brief Computes the maximum-likelihood fit for a selected set of fit parameters.
     *
     * The likelihood backend is selected through @ref StatisticConfig::likelihood_mode.
     * This method also updates the internal state required by contour and scan methods.
     *
     * @param p_specs Ordered list of fit-parameter identifiers.
     * @return Map-based fit result.
     *
     * @throws std::invalid_argument if @p p_specs resolves to an empty fit-parameter set.
     */
    FitResultWithMaps compute_MLE(const std::vector<ParamId>& p_specs);

    /**
     * @brief Computes a two-dimensional confidence contour for the last successful MLE.
     * @param p1 First fit parameter.
     * @param p2 Second fit parameter.
     * @param z Gaussian-equivalent significance level.
     * @param bounds Bounds as {xmin, xmax, ymin, ymax}.
     * @param options Contour and profiling options.
     * @return Extracted contour.
     *
     * @throws std::runtime_error if no MLE has been computed.
     * @throws std::invalid_argument if the parameters are not part of the last fit or are identical.
     */
    Contour confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options);

    /**
     * @brief Prepares a likelihood object for manual scans without running a full MLE.
     *
     * The current central values of the fit parameters and nuisance parameters
     * are stored as the default scan reference point. A later call to
     * compute_MLE(...) or set_manual_scan_point(...) overrides this reference.
     *
     * @param p_specs Ordered list of fit parameters to expose in the scan.
     */
    void prepare_likelihood_for_scan(const std::vector<ParamId>& p_specs);

    /**
     * @brief Sets the reference point used by subsequent likelihood scans.
     * @param p_hat Map of fit-parameter values at the manual reference point.
     * @param eta_hat Map of nuisance values at the manual reference point.
     *
     * @throws std::runtime_error if no likelihood has been prepared.
     * @throws std::invalid_argument if one required parameter value is missing.
     */
    void set_manual_scan_point(const std::map<ParamId, double>& p_hat,
                            const std::map<ParamId, double>& eta_hat);

    /**
     * @brief Evaluates the current likelihood on a regular 2D grid around the active reference point.
     * @param p1 First scanned parameter.
     * @param p2 Second scanned parameter.
     * @param x_half_width Half-width of the scan window along @p p1.
     * @param y_half_width Half-width of the scan window along @p p2.
     * @param nx Number of x grid points.
     * @param ny Number of y grid points.
     * @return Filled scan grid with NLL and delta-NLL values.
     */
    LikelihoodScanGrid scan_likelihood_around_current_point(
        ParamId p1,
        ParamId p2,
        double x_half_width,
        double y_half_width,
        std::size_t nx,
        std::size_t ny
    ) const;

    /**
     * @brief Writes a likelihood scan grid to a CSV file.
     * @param path Output CSV path.
     * @param grid Scan grid to serialize.
     */
    void save_likelihood_scan_csv(const std::string& path,
                                const LikelihoodScanGrid& grid) const;

    /**
     * @brief Updates the full statistical cache for the selected fit parameters.
     * @param p_specs Fit parameters to detach from the nuisance set.
     */
    void update_cache(const std::vector<ParamId>& p_specs = std::vector<ParamId>());

    /** @brief Restricts subsequent statistics to a single experiment name. */
    void select_experiment(const std::string& experiment);

    /** @brief Restricts subsequent statistics to the provided set of experiment names. */
    void select_experiments(const std::set<std::string>& experiments);

    /** @brief Restricts subsequent statistics to the provided vector of experiment names. */
    void select_experiments(const std::vector<std::string>& experiments);

    /** @brief Clears any experiment selection and uses all available experiments. */
    void select_experiments_all();

    /** @return True if an experiment selection is currently active. */
    bool has_experiment_selection() const noexcept;

    /** @return Active experiment selection, or an empty set when all experiments are selected. */
    std::set<std::string> selected_experiments() const;

    /**
     * @brief Restricts subsequent statistics to an explicit list of experimental measurements.
     *
     * This is stricter than @ref select_experiments(): the experiment-name filter keeps every
     * measurement from the selected experiments whose theory observable has been registered,
     * whereas this method keeps only the exact (experiment, observable/bin) entries listed here.
     * It is useful to reproduce SuperIso `myobs.in` / arXiv ancillary observable lists exactly.
     */
    void select_experiment_observables(const std::set<ExperimentObs>& observables);

    /** @brief Vector overload for @ref select_experiment_observables. */
    void select_experiment_observables(const std::vector<ExperimentObs>& observables);

    /** @brief Clears the explicit experimental-observable selection. */
    void select_experiment_observables_all();

    /** @return True if an explicit experimental-observable selection is active. */
    bool has_experiment_observable_selection() const noexcept;

    /** @return Active explicit experimental-observable selection, or an empty set when inactive. */
    std::set<ExperimentObs> selected_experiment_observables() const;

    /**
     * @brief Selects all nuisance dependencies relevant to the current observable set.
     * @return Map of selected nuisance identifiers to central values.
     */
    std::map<ParamId, double> get_all_obss_deps();

    /**
     * @brief Resolves initial fit-parameter values from the parameter proxy.
     * @param p_specs Fit-parameter identifiers.
     * @return Map of parameter identifiers to central values.
     */
    std::map<ParamId, double> get_p_specs(const std::vector<ParamId>& p_specs);

    /** @return Correlation matrix for the currently cached nuisance parameters. */
    std::map<ParamId, std::map<ParamId, double>> get_all_correlations();

    /** @return Correlation matrix for the currently cached experimental observables. */
    std::map<ExperimentObs, std::map<ExperimentObs, double>> get_all_obs_correlations();

    /** @return Experimental central values selected for the current observable/experiment set. */
    std::map<ExperimentObs, double> get_obs_exp();

    /**
     * @brief Prints the current internal cache to standard output for debugging.
     */
    void print_cache();

private:
    /** @brief Rebuilds the merged nuisance registry from default and user registries. */
    void rebuild_merged_nuisance_specs();

    /** @brief Clears cached likelihood, fitter and fit-result state. */
    void invalidate_fit_state();

    /**
     * @brief Finds an effective nuisance specification by parameter id.
     * @param pid Parameter identifier.
     * @return Pointer to the specification, or nullptr if none is known.
     */
    const NuisanceSpec* find_nuisance_spec(const ParamId& pid) const;

    /**
     * @brief Resolves the marginal type for a nuisance parameter, including config overrides.
     * @param pid Nuisance parameter identifier.
     * @return Effective marginal type.
     */
    MarginalType resolve_nuisance_marginal_type(const ParamId& pid) const;

    /**
     * @brief Builds a fit-backend parameter definition for one nuisance parameter.
     * @param pid Nuisance parameter identifier.
     * @param value Central value.
     * @param sigma_hint Nominal uncertainty used as optimizer step hint.
     * @return Fit-backend parameter definition with limits when available.
     */
    fit_app::ParameterDefinition make_nuisance_parameter_definition(const ParamId& pid,
                                                                    double value,
                                                                    double sigma_hint) const;

    /**
     * @brief Builds the marginal-distribution configuration for a nuisance parameter.
     * @param pid Nuisance parameter identifier.
     * @param mt Effective marginal type.
     * @return Marginal configuration compatible with @ref MarginalFactory.
     */                                                                
    MarginalConfig make_nuisance_marginal_config(const ParamId& pid,
                                                 MarginalType mt) const;

    /** @brief Combined experiment-name and explicit-observable filter. */
    bool accepts_experiment_observable(const ExperimentObs& exp_obs) const;

    std::shared_ptr<IModel> obs_int;                       ///< Model interface used for observable predictions.
    std::shared_ptr<IStatCorrelationProxy> pscp;           ///< Proxy providing correlation information.
    std::shared_ptr<IStatParameterProxy> pspp;             ///< Proxy providing parameter values and uncertainties.
    MarginalConfigFactory marginal_config_factory_;        ///< Factory configured with the same parameter port.
    std::shared_ptr<IStatSourcesProxy> sp;                 ///< Proxy resolving source/dependency relationships.
    std::shared_ptr<IStatDependencyPruner> dp;             ///< Dependency pruner used to detach fit parameters.
    std::shared_ptr<INuisanceReader> nuisance_reader_;     ///< Reader for default/user nuisance specifications.
    std::shared_ptr<IStatParamOptimizerProxy> spop;        ///< Optional optimizer proxy retained for workflow integration.

    StatisticConfig config;                                ///< Runtime configuration.
    StatCache cache;                                       ///< Current statistical cache.

    NuisanceRegistry default_nuisance_specs_;              ///< Default nuisance registry.
    NuisanceRegistry user_nuisance_specs_;                 ///< User nuisance registry.
    NuisanceRegistry merged_nuisance_specs_;               ///< Effective nuisance registry after user overrides.
    std::optional<fs::path> current_user_nuisance_file_;   ///< Optional explicit user nuisance file.

    std::shared_ptr<LikelihoodContext> last_ctx_;          ///< Last likelihood context built by MLE or scan preparation.
    std::shared_ptr<BaseLikelihood> last_like_;            ///< Last likelihood object built by MLE or scan preparation.
    std::shared_ptr<MLFitter> last_fitter_;                ///< Last fitter object used for MLE/contours.
    FitResult last_fit_raw_;                               ///< Last raw fit result with vector-based storage.

    std::vector<ParamId> last_fit_param_ids_;              ///< Ordered fit-parameter ids used by the last likelihood.
    std::vector<ParamId> last_nuisance_ids_;               ///< Ordered nuisance ids used by the last likelihood.
    std::map<ParamId, std::size_t> last_fit_param_index_;  ///< Lookup from fit-parameter id to vector index.

    std::vector<double> last_scan_p_;                      ///< Manual or prepared fit-parameter reference values for scans.
    std::vector<double> last_scan_eta_;                    ///< Manual or prepared nuisance reference values for scans.
    bool has_manual_scan_point_ = false;                   ///< True if a manual scan reference point has been provided.

    std::vector<ParamId> last_detached_fit_params_;        ///< Fit parameters detached from the dependency graph during cache update.
    std::vector<std::pair<ParameterType, std::string>> last_detached_fit_blocks_; ///< Fit blocks detached during cache update.

    std::optional<std::set<std::string>> selected_experiments_; ///< Active experiment selection, if any.
    std::optional<std::set<ExperimentObs>> selected_experiment_observables_; ///< Active exact experimental-observable selection, if any.
};

#endif