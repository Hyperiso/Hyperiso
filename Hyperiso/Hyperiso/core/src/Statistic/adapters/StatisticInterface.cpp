#include "StatisticInterface.h"

#include "RuntimeNuisancePathsAdapter.h"

StatisticInterface::StatisticInterface(StatisticConfig config, std::shared_ptr<ObservableInterface> oi_) {
    std::shared_ptr<ObservableInterface> oi = oi_;
    std::shared_ptr<IStatParamOptimizerProxy> spop = std::make_shared<StatParamOptimizerProxy>();
    std::shared_ptr<IModel> oia = std::make_shared<ObservableInterfaceProxy>(oi, spop);
    std::shared_ptr<IStatCorrelationProxy> pscp = std::make_shared<StatCorrelationProxy>();
    std::shared_ptr<IStatParameterProxy> pspp = std::make_shared<StatParameterProxy>();
    std::shared_ptr<IStatSourcesProxy> sp = std::make_shared<StatParamSourcesProxy>();
    std::shared_ptr<IStatDependencyPruner> sdp = std::make_shared<StatDependencyPruner>();
    std::shared_ptr<INuisancePathsProvider> npp =
        std::make_shared<RuntimeNuisancePathsAdapter>();
    std::shared_ptr<INuisanceReader> nr = std::make_shared<NuisanceReader>(npp);

    manager = std::make_shared<StatisticManager>(config, oia, pscp, pspp, sp, sdp, nr, spop);

    manager->select_experiments_all();
    manager->update_cache();
}

void StatisticInterface::select_experiment(const std::string& experiment) {
    manager->select_experiment(experiment);
    manager->update_cache();
}

void StatisticInterface::select_experiments(const std::vector<std::string>& experiments) {
    manager->select_experiments(experiments);
    manager->update_cache();
}

void StatisticInterface::select_experiments_all() {
    manager->select_experiments_all();
    manager->update_cache();
}

bool StatisticInterface::has_experiment_selection() const noexcept {
    return manager->has_experiment_selection();
}

std::set<std::string> StatisticInterface::selected_experiments() const {
    return manager->selected_experiments();
}

void StatisticInterface::select_experiment_observables(
    const std::vector<ExperimentObs>& observables
) {
    manager->select_experiment_observables(observables);
    manager->update_cache();
}

void StatisticInterface::select_experiment_observables_all() {
    manager->select_experiment_observables_all();
    manager->update_cache();
}

bool StatisticInterface::has_experiment_observable_selection() const noexcept {
    return manager->has_experiment_observable_selection();
}

std::set<ExperimentObs> StatisticInterface::selected_experiment_observables() const {
    return manager->selected_experiment_observables();
}

std::map<BinnedObservableId, GaussianSummary> StatisticInterface::compute_uncertainties() {
    return manager->compute_uncertainties();
}

std::map<ParamId, double> StatisticInterface::get_active_observable_dependencies() {
    return manager->get_all_obss_deps();
}

MCResult StatisticInterface::compute_uncertainties_and_sampling() {
    return manager->compute_uncertainties_and_sampling();
}

FitResultWithMaps StatisticInterface::compute_MLE(const std::vector<ParamId>& p_specs) {
    return manager->compute_MLE(p_specs);
}

Contour StatisticInterface::compute_confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options) {
    return manager->confidence_contour(p1, p2, z, bounds, options);
}

void StatisticInterface::reload_nuisance_specs() {
    manager->reload_nuisance_specs();
    manager->update_cache();
}

void StatisticInterface::set_nuisance_user_file(const std::string& user_yaml_path) {
    manager->set_nuisance_user_file(user_yaml_path);
    manager->update_cache();
}

void StatisticInterface::clear_nuisance_user_file() {
    manager->clear_nuisance_user_file();
    manager->update_cache();
}

void StatisticInterface::prepare_likelihood_for_scan(const std::vector<ParamId>& p_specs) {
    manager->prepare_likelihood_for_scan(p_specs);
}

void StatisticInterface::set_manual_scan_point(const std::map<ParamId, double>& p_hat,
                            const std::map<ParamId, double>& eta_hat) {
    manager->set_manual_scan_point(p_hat, eta_hat);
}

LikelihoodScanGrid StatisticInterface::scan_likelihood_around_current_point(
    ParamId p1,
    ParamId p2,
    double x_half_width,
    double y_half_width,
    std::size_t nx,
    std::size_t ny
) const {
    return manager->scan_likelihood_around_current_point(
        p1, p2, x_half_width, y_half_width, nx, ny
    );
}

void StatisticInterface::save_likelihood_scan_csv(const std::string& path,
                                const LikelihoodScanGrid& grid) const {
    manager->save_likelihood_scan_csv(path, grid);
}

void StatisticInterface::update_cache(const std::vector<ParamId>& p_specs) {
    manager->update_cache(p_specs);
}

std::map<ParamId, double>StatisticInterface:: get_all_obss_deps() {
    return manager->get_all_obss_deps();
}

std::map<ParamId, double> StatisticInterface::get_p_specs(const std::vector<ParamId>& p_specs) {
    return manager->get_p_specs(p_specs);
}

std::map<ParamId, std::map<ParamId, double>> StatisticInterface::get_all_correlations() {
    return manager->get_all_correlations();
}

std::map<ExperimentObs, std::map<ExperimentObs, double>> StatisticInterface::get_all_obs_correlations() {
    return manager->get_all_obs_correlations();
}

std::map<ExperimentObs, double> StatisticInterface::get_obs_exp() {
    return manager->get_obs_exp();
}

void StatisticInterface::print_cache() {
    manager->print_cache();
}
