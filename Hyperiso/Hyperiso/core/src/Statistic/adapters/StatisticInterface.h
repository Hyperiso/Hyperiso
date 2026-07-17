#ifndef STATISTIC_INTERFACE_H
#define STATISTIC_INTERFACE_H

#include "StatisticManager.h"
#include "ObservableInterfaceProxy.h"
#include "StatCorrelationProxy.h"
#include "StatParameterProxy.h"
#include "ObservableInterface.h"
#include "StatParamSourcesProxy.h"
#include "StatDependencyPruner.h"
#include "NuisanceReader.h"

/**
 * @file StatisticInterface.h
 * @brief High-level facade for statistical computations.
 *
 * This interface wires the observable layer, parameter/correlation proxies,
 * nuisance readers and @ref StatisticManager into a single user-facing class.
 */

/**
 * @class StatisticInterface
 * @brief Facade exposing statistical workflows for an ObservableInterface.
 *
 * The class owns a configured @ref StatisticManager and forwards high-level
 * operations such as uncertainty propagation, maximum-likelihood fits,
 * confidence contours and likelihood scans.
 */
class StatisticInterface {
public:
    /**
     * @brief Constructs the statistical interface around an observable backend.
     *
     * The constructor creates the default proxy stack, initializes the
     * underlying @ref StatisticManager, selects all experiments and builds the
     * initial cache.
     *
     * @param config Statistic configuration.
     * @param oi_    Observable interface to use for predictions.
     */
    StatisticInterface(StatisticConfig config, std::shared_ptr<ObservableInterface> oi_);

    /**
     * @brief Restricts the analysis to a single experiment.
     *
     * @param experiment Experiment name to select.
     */
    void select_experiment(const std::string& experiment);

    /**
     * @brief Restricts the analysis to a set of experiments.
     *
     * @param experiments Experiment names to select.
     */
    void select_experiments(const std::vector<std::string>& experiments);

    /**
     * @brief Clears any experiment restriction and uses all experiments.
     */
    void select_experiments_all();

    /**
     * @brief Tests whether an experiment selection is active.
     *
     * @return True if only a subset of experiments is selected.
     */
    bool has_experiment_selection() const noexcept;

    /**
     * @brief Returns the currently selected experiments.
     *
     * @return Selected experiment names, or an empty set when all experiments
     *         are selected.
     */
    std::set<std::string> selected_experiments() const;

    /**
     * @brief Restricts the analysis to exact experiment/observable/bin entries.
     *
     * @param observables Exact experimental measurements to keep.
     */
    void select_experiment_observables(const std::vector<ExperimentObs>& observables);

    /**
     * @brief Clears any exact experimental-observable restriction.
     */
    void select_experiment_observables_all();

    /**
     * @brief Tests whether an exact experimental-observable selection is active.
     *
     * @return True if exact measurements were selected explicitly.
     */
    bool has_experiment_observable_selection() const noexcept;

    /**
     * @brief Returns the active exact experimental-observable selection.
     *
     * @return Selected measurements, or an empty set when no exact selection is active.
     */
    std::set<ExperimentObs> selected_experiment_observables() const;

    /**
     * @brief Computes Gaussian uncertainty summaries for all active observables.
     *
     * @return Map from binned observable id to Gaussian summary.
     */
    std::map<BinnedObservableId, GaussianSummary> compute_uncertainties();

    /**
     * @brief Returns the nuisance/input parameters currently seen by the statistical manager.
     *
     * This is mostly useful for dev/debug workflows. It runs the same observable
     * dependency path as the cache builder: IModel::get_obs_deps(), leaf-source
     * expansion through the source proxy, and the manager's sensitivity filtering.
     * Runtime/lambda observables are therefore included as long as their
     * dependencies were declared on the ObservableInterface side.
     */
    std::map<ParamId, double> get_active_observable_dependencies();

    /**
     * @brief Runs uncertainty propagation and returns both samples and summaries.
     *
     * @return Monte Carlo realization, Gaussian summaries and covariance data.
     */
    MCResult compute_uncertainties_and_sampling();
    
    /**
     * @brief Runs a maximum-likelihood fit for the requested fit parameters.
     *
     * @param p_specs Fit parameters to float.
     *
     * @return Fit result keyed by parameter identifiers.
     */
    FitResultWithMaps compute_MLE(const std::vector<ParamId>& p_specs);

    /**
     * @brief Computes a two-dimensional confidence contour.
     *
     * @param p1      First fit parameter.
     * @param p2      Second fit parameter.
     * @param z       Gaussian-equivalent confidence level in standard deviations.
     * @param bounds  Contour bounds as {xmin, xmax, ymin, ymax}.
     * @param options Contour and profiling options.
     *
     * @return Extracted contour.
     */
    Contour compute_confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options);

    /**
     * @brief Reloads default and user nuisance specifications.
     */
    void reload_nuisance_specs();

    /**
     * @brief Uses a specific user nuisance configuration file.
     *
     * @param user_yaml_path Path to the user YAML file.
     */
    void set_nuisance_user_file(const std::string& user_yaml_path);

    /**
     * @brief Restores the default user nuisance file lookup.
     */
    void clear_nuisance_user_file();

    /**
     * @brief Builds and stores the likelihood objects required for scan calls.
     *
     * The current central values are also stored as the default scan center,
     * so lightweight scans can run without computing a full MLE first.
     * Calling compute_MLE(...) or set_manual_scan_point(...) later replaces
     * that reference point.
     *
     * @param p_specs Fit parameters defining the scan parameter space.
     */
    void prepare_likelihood_for_scan(const std::vector<ParamId>& p_specs);

    /**
     * @brief Manually sets the central point used by likelihood scans.
     *
     * @param p_hat   Fit-parameter values at the scan center.
     * @param eta_hat Nuisance-parameter values at the scan center.
     */
    void set_manual_scan_point(const std::map<ParamId, double>& p_hat,
                               const std::map<ParamId, double>& eta_hat);

    /**
     * @brief Evaluates a rectangular likelihood scan around the current point.
     *
     * @param p1           First scan parameter.
     * @param p2           Second scan parameter.
     * @param x_half_width Half-width along the first parameter axis.
     * @param y_half_width Half-width along the second parameter axis.
     * @param nx           Number of grid points along x.
     * @param ny           Number of grid points along y.
     *
     * @return Likelihood scan grid.
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
     * @brief Saves a likelihood scan grid as a CSV file.
     *
     * @param path Output CSV path.
     * @param grid Grid to serialize.
     */
    void save_likelihood_scan_csv(const std::string& path,
                                  const LikelihoodScanGrid& grid) const;

    /**
     * @brief Refreshes the cached statistical inputs.
     *
     * @param p_specs Optional fit-parameter list used to define nuisance
     *                pruning and dependency handling.
     */
    void update_cache(const std::vector<ParamId>& p_specs = std::vector<ParamId>());

    /**
     * @copydoc StatisticManager::get_all_obss_deps()
     */
    std::map<ParamId, double> get_all_obss_deps();

    /**
     * @copydoc StatisticManager::get_p_specs(const std::vector<ParamId>&)
     */
    std::map<ParamId, double> get_p_specs(const std::vector<ParamId>& p_specs);

    /**
     * @copydoc StatisticManager::get_all_correlations()
     */
    std::map<ParamId, std::map<ParamId, double>> get_all_correlations();

    /**
     * @copydoc StatisticManager::get_all_obs_correlations()
     */
    std::map<ExperimentObs, std::map<ExperimentObs, double>> get_all_obs_correlations();

    /**
     * @copydoc StatisticManager::get_obs_exp()
     */
    std::map<ExperimentObs, double> get_obs_exp();

    /**
     * @brief Prints the current manager cache for diagnostics.
     */
    void print_cache();

private:
    std::shared_ptr<StatisticManager> manager;  ///< Underlying statistical workflow manager.
};

#endif