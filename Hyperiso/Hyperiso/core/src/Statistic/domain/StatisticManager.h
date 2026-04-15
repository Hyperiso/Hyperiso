// #ifndef __STATISTIC_MANAGER_H__
// #define __STATISTIC_MANAGER_H__

// #include <vector>
// #include <iomanip>
// #include <iostream>
// #include <filesystem>
// #include <optional>

// #include "Include.h"
// #include "IModel.h"
// #include "IStatCorrelationProxy.h"
// #include "IStatParameterProxy.h"
// #include "IStatSourcesProxy.h"
// #include "IStatDependencyPruner.h"
// #include "CovarianceTransformer.h"
// #include "JointDistribution.h"
// #include "RvgNuisanceSampler.h"
// #include "MCEngine.h"
// #include "DistributionFactory.h"
// #include "CopulaFactory.h"
// #include "Fit.h"
// #include "MarginalConfigFactory.h"
// #include "INuisanceReader.h"
// #include "NuisanceSpec.h"


// struct StatisticConfig {
//     std::map<ParamId, MarginalType> override_nuisance_marginals {};
//     std::map<ExperimentObs, MarginalType> override_exp_data_marginals {};
//     CopulaType nuisance_copula_type = CopulaType::GAUSSIAN;
//     CopulaType exp_data_copula_type = CopulaType::GAUSSIAN;
//     std::size_t MC_draws = 100;
//     double skew_abs_threshold=0.2;

//     std::size_t MLE_max_iter = 500;
//     double MLE_tol = 1e-6;

//     double nuisance_relevance_cutoff = 1e-2;
// };

// struct FitResultWithMaps {
//     bool fit_ok {false};
//     std::map<ParamId, double> p_hat;
//     std::map<ParamId, double> eta_hat;
//     std::map<ParamId, double> p_hat_std;
//     std::map<ParamId, std::map<ParamId, double>> p_correlations;
//     double ell_hat {0.0};
// };

// struct StatCache {
//     std::map<ExperimentObs, double> exp_obs;
//     std::map<ParamId, double> eta_specs_real;
//     std::map<ParamId, double> p_specs;
//     std::map<ParamId, std::map<ParamId, double>> SigmaEta;
//     std::map<ExperimentObs, std::map<ExperimentObs, double>> SigmaObs;

//     FitResultWithMaps mle_result;
// };

// class StatisticManager {
// public:
//     StatisticManager(StatisticConfig config,
//                      std::shared_ptr<IModel> obs_int,
//                      std::shared_ptr<IStatCorrelationProxy> pscp,
//                      std::shared_ptr<IStatParameterProxy> pspp,
//                      std::shared_ptr<IStatSourcesProxy> sp,
//                      std::shared_ptr<IStatDependencyPruner> dp,
//                      std::shared_ptr<INuisanceReader> nuisance_reader);

//     std::vector<std::unique_ptr<IMarginalDistribution>> build_nuisance_marginal_distributions();
//     std::unique_ptr<JointDistribution> build_nuisance_distribution();
//     std::unique_ptr<JointDistribution> build_exp_data_distribution();

//     void reload_nuisance_specs();
//     void set_nuisance_user_file(const fs::path& user_yaml_path);
//     void clear_nuisance_user_file();

//     const NuisanceRegistry& default_nuisance_specs() const { return default_nuisance_specs_; }
//     const NuisanceRegistry& user_nuisance_specs() const { return user_nuisance_specs_; }
//     const NuisanceRegistry& merged_nuisance_specs() const { return merged_nuisance_specs_; }


//     std::map<BinnedObservableId, GaussianSummary> compute_uncertainties() {
//         auto sums = this->compute_uncertainties_and_sampling();

//         std::ofstream fs;
//         fs.open("samples.csv");        

//         std::map<BinnedObservableId, GaussianSummary> out;
//         for (const auto& gs : sums.summary) {
//             out[gs.id] = gs;
//             std::cout << gs << std::endl;
//         }
//         return out;
//     }

//     MCResult compute_uncertainties_and_sampling() {
//         update_cache();
//         auto rvg = build_nuisance_distribution();
//         std::vector<ParamId> nuisance_ids = unzip(cache.eta_specs_real).ids;
//         RvgNuisanceSampler sampler(nuisance_ids, std::move(rvg));
//         MonteCarloEngine mc(this->obs_int, sampler, {this->config.MC_draws, this->config.skew_abs_threshold});

//         auto sums = mc.summarize(this->cache.p_specs);

//         return sums;
//     }
    
//     FitResultWithMaps compute_MLE(const std::vector<ParamId>& p_specs);

//     Contour confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options);

//     void update_cache(const std::vector<ParamId>& p_specs = std::vector<ParamId>()) {
//         cache.p_specs = this->get_p_specs(p_specs);
//         for (const auto& [pid, _] : cache.p_specs)
//             dp->detach_parameter(pid.type.value(), pid.block, pid.code);
//         cache.eta_specs_real = this->get_all_obss_deps();
//         for (const auto& [pid, _] : cache.p_specs)
//             cache.eta_specs_real.erase(pid);
//         cache.SigmaEta = this->get_all_correlations();
//         cache.exp_obs = this->get_obs_exp();
//         cache.SigmaObs = this->get_all_obs_correlations();
//     }

//     // std::map<ParamId, double> get_all_obss_deps() {
//     //     std::unordered_set<ParamId> eta_infos;

//     //     for (const auto& obsId : obs_int->get_obs_ids()) {
//     //         for (auto paramId : obs_int->get_obs_deps(obsId.s)) {
//     //             eta_infos.insert(paramId);
//     //         }
//     //     }

//     //     std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);
//     //     std::map<ParamId, double> eta_specs_real_leaf;
//     //     std::map<ParamId, double> delta_rel;

//     //     for (auto& paramId : eta_infos_leaf) {
//     //         double u = pspp->get_param(paramId)->get_combined_std().real();
//     //         double val = pspp->get_param(paramId)->get_val().real();
//     //         if (fpeq(val, 0.0) && !fpeq(u, 0.0)) {
//     //             eta_specs_real_leaf[paramId] = val;
//     //             continue;
//     //         }
                
//     //         if (!std::isfinite(u) || fpeq(std::abs(u), 0.0))
//     //             delta_rel[paramId] = 0;
//     //         else
//     //             delta_rel[paramId] = std::abs(u / val);
//     //     }
            
        
//     //     using T = std::pair<ParamId, double>;
//     //     double delta_rel_max = std::max_element(
//     //         delta_rel.begin(), delta_rel.end(), 
//     //         [] (const T& p, const T& q) { return p.second < q.second; }
//     //     )->second;

//     //     for (auto& paramId : eta_infos_leaf) {
//     //         LOG_INFO("Compared to max relative uncertainty", paramId, delta_rel[paramId] / delta_rel_max);
//     //     }

//     //     for (auto& paramId : eta_infos_leaf) {
//     //         if (delta_rel[paramId] / delta_rel_max > config.nuisance_relevance_cutoff)
//     //             eta_specs_real_leaf[paramId] = pspp->get_param(paramId)->get_val();
//     //     }

//     //     LOG_INFO("Significant nuisances");
//     //     for (const auto& [pid, val] : eta_specs_real_leaf) {
//     //         LOG_INFO(pid, val);
//     //     }

//     //     return eta_specs_real_leaf;
//     // }

//     // std::map<ParamId, double> get_p_specs(const std::vector<ParamId>& p_specs) {
//     //     std::map<ParamId, double> out;
//     //     for (auto elem : p_specs) {            
//     //         out[elem] = pspp->get_param(elem)->get_val();
//     //     }
//     //     return out;
//     // }
//     // std::map<ParamId, std::map<ParamId, double>> get_all_correlations() {
//     //     std::map<ParamId, std::map<ParamId, double>> res;
//     //     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);
//     //     res = ct.transform(this->cache.eta_specs_real);
//     //     return res;
//     // }
    
//     // std::map<ExperimentObs, std::map<ExperimentObs, double>> get_all_obs_correlations() {
//     //     std::map<ExperimentObs, std::map<ExperimentObs, double>> res;
//     //     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);

//     //     res = ct.transform(this->cache.exp_obs);
//     //     return res;
//     // }

//     // std::map<ExperimentObs, double> get_obs_exp() {
//     //     std::map<ExperimentObs, double> out;

//     //     for (const auto& obsId : obs_int->get_obs_ids()) {
//     //         auto _ = pspp->get_obs_param(obsId);
//     //         for (auto _2 : _) {
//     //             out[_2.first] = _2.second->get_val(); 
//     //         }
//     //         // out[obsId] = pspp->get_obs_param(obsId)->get_val();
//     //     }
//     //     return out;
//     // }

//     std::map<ParamId, double> get_all_obss_deps();
//     std::map<ParamId, double> get_p_specs(const std::vector<ParamId>& p_specs);
//     std::map<ParamId, std::map<ParamId, double>> get_all_correlations();
//     std::map<ExperimentObs, std::map<ExperimentObs, double>> get_all_obs_correlations();
//     std::map<ExperimentObs, double> get_obs_exp();

//     void print_cache();

// private:
//     void rebuild_merged_nuisance_specs();
//     void invalidate_fit_state();

//     const NuisanceSpec* find_nuisance_spec(const ParamId& pid) const;
//     MarginalType resolve_nuisance_marginal_type(const ParamId& pid) const;

//     fit_app::ParameterDefinition make_nuisance_parameter_definition(const ParamId& pid,
//                                                                     double value,
//                                                                     double sigma_hint) const;

//     MarginalConfig make_nuisance_marginal_config(const ParamId& pid,
//                                                  MarginalType mt) const;

//     std::shared_ptr<IModel> obs_int;
//     std::shared_ptr<IStatCorrelationProxy> pscp;
//     std::shared_ptr<IStatParameterProxy> pspp;
//     std::shared_ptr<IStatSourcesProxy> sp;
//     std::shared_ptr<IStatDependencyPruner> dp;
//     std::shared_ptr<INuisanceReader> nuisance_reader_;

//     StatisticConfig config;
//     StatCache cache;
    
//     NuisanceRegistry default_nuisance_specs_;
//     NuisanceRegistry user_nuisance_specs_;
//     NuisanceRegistry merged_nuisance_specs_;
//     std::optional<fs::path> current_user_nuisance_file_;

//     std::shared_ptr<LikelihoodContext> last_ctx_;
//     std::shared_ptr<BaseLikelihood> last_like_;
//     std::shared_ptr<MLFitter> last_fitter_;
//     FitResult last_fit_raw_;

//     std::vector<ParamId> last_fit_param_ids_;
//     std::vector<ParamId> last_nuisance_ids_;
//     std::map<ParamId, std::size_t> last_fit_param_index_;
// };

// #endif

#ifndef __STATISTIC_MANAGER_H__
#define __STATISTIC_MANAGER_H__

#include <vector>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <optional>

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
#include "DistributionFactory.h"
#include "CopulaFactory.h"
#include "Fit.h"
#include "MarginalConfigFactory.h"
#include "INuisanceReader.h"
#include "NuisanceSpec.h"


// struct StatisticConfig {
//     std::map<ParamId, MarginalType> override_nuisance_marginals {};
//     std::map<ExperimentObs, MarginalType> override_exp_data_marginals {};
//     CopulaType nuisance_copula_type = CopulaType::GAUSSIAN;
//     CopulaType exp_data_copula_type = CopulaType::GAUSSIAN;
//     std::size_t MC_draws = 100;
//     double skew_abs_threshold=0.2;

//     std::size_t MLE_max_iter = 500;
//     double MLE_tol = 1e-6;
//     unsigned MLE_strategy = 2;
//     bool MLE_run_hesse = true;
//     bool MLE_request_minos = false;
//     bool MLE_verbose = false;

//     double nuisance_relevance_cutoff = 1e-2;
// };

struct StatisticConfig {
    std::map<ParamId, MarginalType> override_nuisance_marginals {};
    std::map<ExperimentObs, MarginalType> override_exp_data_marginals {};
    CopulaType nuisance_copula_type = CopulaType::GAUSSIAN;
    CopulaType exp_data_copula_type = CopulaType::GAUSSIAN;
    std::size_t MC_draws = 100;
    double skew_abs_threshold = 0.2;

    std::size_t MLE_max_iter = 500;
    double MLE_tol = 1e-6;
    unsigned MLE_strategy = 2;
    bool MLE_run_hesse = true;
    bool MLE_request_minos = false;
    bool MLE_verbose = false;

    double nuisance_relevance_cutoff = 1e-2;

    // Nouveau : pruning par sensibilité locale du modèle
    bool nuisance_sensitivity_pruning = true;

    // On teste eta_j -> eta_j ± probe_sigmas * sigma_j
    double nuisance_sensitivity_probe_sigmas = 1.0;

    // Une nuisance est gardée si son effet 1σ dépasse au moins
    // un des deux seuils ci-dessous.
    double nuisance_sensitivity_rel_cutoff = 1e-4;
    double nuisance_sensitivity_abs_cutoff = 1e-10;

    // Evite de diviser par des échelles trop petites quand une observable ~ 0
    double nuisance_sensitivity_scale_floor = 1e-3;

    bool MLE_trace_first_evals = false;
    std::size_t MLE_trace_max_evals = 25;

    bool MLE_allow_profile_hessian_fallback = true;
    double MLE_profile_hessian_step_scale = 1.0;
    double MLE_profile_hessian_eig_floor_rel = 1e-8;

};

struct LikelihoodScanPoint {
    double x = 0.0;
    double y = 0.0;
    double nll = 0.0;
    double delta_nll = 0.0;
};

struct LikelihoodScanGrid {
    ParamId x_param;
    ParamId y_param;
    double x_center = 0.0;
    double y_center = 0.0;
    std::size_t nx = 0;
    std::size_t ny = 0;
    std::vector<LikelihoodScanPoint> points;
};

struct FitResultWithMaps {
    bool fit_ok {false};
    std::map<ParamId, double> p_hat;
    std::map<ParamId, double> eta_hat;
    std::map<ParamId, double> p_hat_std;
    std::map<ParamId, std::map<ParamId, double>> p_correlations;
    double ell_hat {0.0};
};

struct StatCache {
    std::map<ExperimentObs, double> exp_obs;
    std::map<ParamId, double> eta_specs_real;
    std::map<ParamId, double> p_specs;
    std::map<ParamId, std::map<ParamId, double>> SigmaEta;
    std::map<ExperimentObs, std::map<ExperimentObs, double>> SigmaObs;

    FitResultWithMaps mle_result;
};

class StatisticManager {
public:
    StatisticManager(StatisticConfig config,
                     std::shared_ptr<IModel> obs_int,
                     std::shared_ptr<IStatCorrelationProxy> pscp,
                     std::shared_ptr<IStatParameterProxy> pspp,
                     std::shared_ptr<IStatSourcesProxy> sp,
                     std::shared_ptr<IStatDependencyPruner> dp,
                     std::shared_ptr<INuisanceReader> nuisance_reader);

    std::vector<std::unique_ptr<IMarginalDistribution>> build_nuisance_marginal_distributions();
    std::unique_ptr<JointDistribution> build_nuisance_distribution();
    std::unique_ptr<JointDistribution> build_exp_data_distribution();

    void reload_nuisance_specs();
    void set_nuisance_user_file(const fs::path& user_yaml_path);
    void clear_nuisance_user_file();

    const NuisanceRegistry& default_nuisance_specs() const { return default_nuisance_specs_; }
    const NuisanceRegistry& user_nuisance_specs() const { return user_nuisance_specs_; }
    const NuisanceRegistry& merged_nuisance_specs() const { return merged_nuisance_specs_; }


    std::map<BinnedObservableId, GaussianSummary> compute_uncertainties() {
        auto sums = this->compute_uncertainties_and_sampling();

        for (auto elem : this->merged_nuisance_specs_) {
            std::cout << elem.first << std::endl;
            std::cout << elem.second << std::endl;
            std::cout << "-------------------" << std::endl;
        }
        std::ofstream fs;
        fs.open("samples.csv");        

        std::map<BinnedObservableId, GaussianSummary> out;
        for (const auto& gs : sums.summary) {
            out[gs.id] = gs;
            std::cout << gs << std::endl;
        }
        return out;
    }

    MCResult compute_uncertainties_and_sampling() {
        update_cache();
        auto rvg = build_nuisance_distribution();
        std::vector<ParamId> nuisance_ids = unzip(cache.eta_specs_real).ids;
        RvgNuisanceSampler sampler(nuisance_ids, std::move(rvg));
        MonteCarloEngine mc(this->obs_int, sampler, {this->config.MC_draws, this->config.skew_abs_threshold});

        auto sums = mc.summarize(this->cache.p_specs);

        return sums;
    }
    
    FitResultWithMaps compute_MLE(const std::vector<ParamId>& p_specs);

    Contour confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, ContourOptions options);

    void prepare_likelihood_for_scan(const std::vector<ParamId>& p_specs);

    void set_manual_scan_point(const std::map<ParamId, double>& p_hat,
                            const std::map<ParamId, double>& eta_hat);

    LikelihoodScanGrid scan_likelihood_around_current_point(
        ParamId p1,
        ParamId p2,
        double x_half_width,
        double y_half_width,
        std::size_t nx,
        std::size_t ny
    ) const;

    void save_likelihood_scan_csv(const std::string& path,
                                const LikelihoodScanGrid& grid) const;

    // void update_cache(const std::vector<ParamId>& p_specs = std::vector<ParamId>()) {
    //     cache.p_specs = this->get_p_specs(p_specs);
    //     for (const auto& [pid, _] : cache.p_specs) {
    //     if (!pid.type.has_value()) {
    //         continue;
    //     }

    //     dp->detach_block(pid.type.value(), pid.block);

    //     dp->detach_parameter(pid.type.value(), pid.block, pid.code);
    // }

    //     cache.eta_specs_real = this->get_all_obss_deps();
    //     for (const auto& [pid, _] : cache.p_specs)
    //         cache.eta_specs_real.erase(pid);
    //     cache.SigmaEta = this->get_all_correlations();
    //     cache.exp_obs = this->get_obs_exp();
    //     cache.SigmaObs = this->get_all_obs_correlations();
    // }

    void update_cache(const std::vector<ParamId>& p_specs = std::vector<ParamId>()) {
        for (const auto& [tp, block] : last_detached_fit_blocks_) {
            dp->reattach_block(tp, block);
        }
        for (const auto& pid : last_detached_fit_params_) {
            if (pid.type.has_value()) {
                dp->reattach_parameter(pid.type.value(), pid.block, pid.code);
            }
        }
        last_detached_fit_blocks_.clear();
        last_detached_fit_params_.clear();

        cache.p_specs = this->get_p_specs(p_specs);

        std::unordered_set<std::string> seen_blocks;

        for (const auto& [pid, _] : cache.p_specs) {
            if (!pid.type.has_value()) {
                continue;
            }

            const auto tp = pid.type.value();
            const std::string block_key =
                std::to_string(static_cast<int>(tp)) + "::" + pid.block.to_string();

            if (!seen_blocks.contains(block_key)) {
                dp->detach_block(tp, pid.block);
                last_detached_fit_blocks_.push_back({tp, pid.block});
                seen_blocks.insert(block_key);
            }

            dp->detach_parameter(tp, pid.block, pid.code);
            last_detached_fit_params_.push_back(pid);
        }

        cache.eta_specs_real = this->get_all_obss_deps();
        for (const auto& [pid, _] : cache.p_specs)
            cache.eta_specs_real.erase(pid);

        cache.SigmaEta = this->get_all_correlations();
        cache.exp_obs = this->get_obs_exp();
        cache.SigmaObs = this->get_all_obs_correlations();
    }

    // std::map<ParamId, double> get_all_obss_deps() {
    //     std::unordered_set<ParamId> eta_infos;

    //     for (const auto& obsId : obs_int->get_obs_ids()) {
    //         for (auto paramId : obs_int->get_obs_deps(obsId.s)) {
    //             eta_infos.insert(paramId);
    //         }
    //     }

    //     std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);
    //     std::map<ParamId, double> eta_specs_real_leaf;
    //     std::map<ParamId, double> delta_rel;

    //     for (auto& paramId : eta_infos_leaf) {
    //         double u = pspp->get_param(paramId)->get_combined_std().real();
    //         double val = pspp->get_param(paramId)->get_val().real();
    //         if (fpeq(val, 0.0) && !fpeq(u, 0.0)) {
    //             eta_specs_real_leaf[paramId] = val;
    //             continue;
    //         }
                
    //         if (!std::isfinite(u) || fpeq(std::abs(u), 0.0))
    //             delta_rel[paramId] = 0;
    //         else
    //             delta_rel[paramId] = std::abs(u / val);
    //     }
            
        
    //     using T = std::pair<ParamId, double>;
    //     double delta_rel_max = std::max_element(
    //         delta_rel.begin(), delta_rel.end(), 
    //         [] (const T& p, const T& q) { return p.second < q.second; }
    //     )->second;

    //     for (auto& paramId : eta_infos_leaf) {
    //         LOG_INFO("Compared to max relative uncertainty", paramId, delta_rel[paramId] / delta_rel_max);
    //     }

    //     for (auto& paramId : eta_infos_leaf) {
    //         if (delta_rel[paramId] / delta_rel_max > config.nuisance_relevance_cutoff)
    //             eta_specs_real_leaf[paramId] = pspp->get_param(paramId)->get_val();
    //     }

    //     LOG_INFO("Significant nuisances");
    //     for (const auto& [pid, val] : eta_specs_real_leaf) {
    //         LOG_INFO(pid, val);
    //     }

    //     return eta_specs_real_leaf;
    // }

    // std::map<ParamId, double> get_p_specs(const std::vector<ParamId>& p_specs) {
    //     std::map<ParamId, double> out;
    //     for (auto elem : p_specs) {            
    //         out[elem] = pspp->get_param(elem)->get_val();
    //     }
    //     return out;
    // }
    // std::map<ParamId, std::map<ParamId, double>> get_all_correlations() {
    //     std::map<ParamId, std::map<ParamId, double>> res;
    //     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);
    //     res = ct.transform(this->cache.eta_specs_real);
    //     return res;
    // }
    
    // std::map<ExperimentObs, std::map<ExperimentObs, double>> get_all_obs_correlations() {
    //     std::map<ExperimentObs, std::map<ExperimentObs, double>> res;
    //     CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);

    //     res = ct.transform(this->cache.exp_obs);
    //     return res;
    // }

    // std::map<ExperimentObs, double> get_obs_exp() {
    //     std::map<ExperimentObs, double> out;

    //     for (const auto& obsId : obs_int->get_obs_ids()) {
    //         auto _ = pspp->get_obs_param(obsId);
    //         for (auto _2 : _) {
    //             out[_2.first] = _2.second->get_val(); 
    //         }
    //         // out[obsId] = pspp->get_obs_param(obsId)->get_val();
    //     }
    //     return out;
    // }

    std::map<ParamId, double> get_all_obss_deps();
    std::map<ParamId, double> get_p_specs(const std::vector<ParamId>& p_specs);
    std::map<ParamId, std::map<ParamId, double>> get_all_correlations();
    std::map<ExperimentObs, std::map<ExperimentObs, double>> get_all_obs_correlations();
    std::map<ExperimentObs, double> get_obs_exp();

    void print_cache();

private:
    void rebuild_merged_nuisance_specs();
    void invalidate_fit_state();

    const NuisanceSpec* find_nuisance_spec(const ParamId& pid) const;
    MarginalType resolve_nuisance_marginal_type(const ParamId& pid) const;

    fit_app::ParameterDefinition make_nuisance_parameter_definition(const ParamId& pid,
                                                                    double value,
                                                                    double sigma_hint) const;

    MarginalConfig make_nuisance_marginal_config(const ParamId& pid,
                                                 MarginalType mt) const;

    std::shared_ptr<IModel> obs_int;
    std::shared_ptr<IStatCorrelationProxy> pscp;
    std::shared_ptr<IStatParameterProxy> pspp;
    std::shared_ptr<IStatSourcesProxy> sp;
    std::shared_ptr<IStatDependencyPruner> dp;
    std::shared_ptr<INuisanceReader> nuisance_reader_;

    StatisticConfig config;
    StatCache cache;
    
    NuisanceRegistry default_nuisance_specs_;
    NuisanceRegistry user_nuisance_specs_;
    NuisanceRegistry merged_nuisance_specs_;
    std::optional<fs::path> current_user_nuisance_file_;

    std::shared_ptr<LikelihoodContext> last_ctx_;
    std::shared_ptr<BaseLikelihood> last_like_;
    std::shared_ptr<MLFitter> last_fitter_;
    FitResult last_fit_raw_;

    std::vector<ParamId> last_fit_param_ids_;
    std::vector<ParamId> last_nuisance_ids_;
    std::map<ParamId, std::size_t> last_fit_param_index_;

    std::vector<double> last_scan_p_;
    std::vector<double> last_scan_eta_;
    bool has_manual_scan_point_ = false;

    std::vector<ParamId> last_detached_fit_params_;
    std::vector<std::pair<ParameterType, std::string>> last_detached_fit_blocks_;
};

#endif