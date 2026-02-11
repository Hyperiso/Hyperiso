#ifndef __STATISTIC_MANAGER_H__
#define __STATISTIC_MANAGER_H__

#include <vector>
#include <iomanip>
#include <iostream>
#include "Include.h"
#include "IModel.h"
#include "IStatCorrelationProxy.h"
#include "IStatParameterProxy.h"
#include "IStatSourcesProxy.h"
#include "CovarianceTransformer.h"
#include "JointDistribution.h"
#include "RvgNuisanceSampler.h"
#include "MCEngine.h"
#include "MarginalFactory.h"
#include "CopulaFactory.h"
#include "Fit.h"
#include "MarginalConfigFactory.h"

struct StatisticConfig {
    std::map<ObservableId, QCDOrder> obss;
    std::vector<ParamId> p_specs;

    std::map<ParamId, MarginalType> override_nuisance_marginals {};
    std::map<ObservableId, MarginalType> override_exp_data_marginals {};
    CopulaType nuisance_copula_type = CopulaType::GAUSSIAN;
    CopulaType exp_data_copula_type = CopulaType::GAUSSIAN;
    std::size_t MC_draws = 100;
    double skew_abs_threshold=0.2;

    std::size_t MLE_max_iter = 500;
    double MLE_tol = 1e-6;
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
    std::map<ObservableId, double> exp_obs;
    std::map<ParamId, double> eta_specs_real;
    std::map<ParamId, double> p_specs;
    std::map<ParamId, std::map<ParamId, double>> SigmaEta;
    std::map<ObservableId, std::map<ObservableId, double>> SigmaObs;

    FitResultWithMaps mle_result;
};

enum class CLMethod {
    SLICE,
    PROJECT
};

class StatisticManager {
public:
    StatisticManager(StatisticConfig config, std::shared_ptr<IModel> obs_int, 
        std::shared_ptr<IStatCorrelationProxy> pscp, std::shared_ptr<IStatParameterProxy> pspp,
        std::shared_ptr<IStatSourcesProxy> sp) 
    : config(config), obs_int(obs_int), pscp(pscp), pspp(pspp), sp(sp) {
        obs_int->add_observables(config.obss);
    }

    std::vector<std::unique_ptr<IMarginalDistribution>> build_nuisance_marginal_distributions();
    std::unique_ptr<JointDistribution> build_nuisance_distribution();
    std::unique_ptr<JointDistribution> build_exp_data_distribution();

    std::map<ObservableId, GaussianSummary> compute_uncertainties() {
        auto sums = this->compute_uncertainties_and_sampling();

        // Debug print
        for (auto sum : sums.summary) {
            std::cout << sum << std::endl;
        }

        return zip(unzip(config.obss).ids, sums.summary);
    }

    MCResult compute_uncertainties_and_sampling() {
        auto rvg = build_nuisance_distribution();
        std::vector<ParamId> nuisance_ids = unzip(cache.eta_specs_real).ids;
        RvgNuisanceSampler sampler(nuisance_ids, std::move(rvg));
        MonteCarloEngine mc(this->obs_int, sampler, {this->config.MC_draws, this->config.skew_abs_threshold});

        auto sums = mc.summarize(this->cache.p_specs);

        return sums;
    }
    
    FitResultWithMaps compute_MLE() {

        // Build Likelihood context

        auto unzipped_fit_params = unzip(cache.p_specs);
        auto unzipped_nuisances = unzip(cache.eta_specs_real);
        auto unzipped_exp_obs = unzip(cache.exp_obs);

        std::vector<ParamId> p_ids = unzipped_fit_params.ids;
        std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
        std::vector<ObservableId> obs_ids = unzipped_exp_obs.ids;
        
        LikelihoodContext ctx;
        ctx.nuisance_dist = std::move(build_nuisance_distribution());
        ctx.exp_obs_dist = std::move(build_exp_data_distribution());
        ctx.nuisance_central_values = unzipped_nuisances.vals;
        ctx.exp_obs_values = unzipped_exp_obs.vals;

        auto model_fn = [this, obs_ids, p_ids, eta_ids] (const Vec& p_vec, const Vec& eta_vec) -> Vec {
            auto pred_map = this->obs_int->predict_optimized(zip(p_ids, p_vec), zip(eta_ids, eta_vec));
            return unzip(pred_map).vals;
        };

        MLEstimator est(std::move(ctx), model_fn, this->config.MLE_max_iter, this->config.MLE_tol);

        FitResult fr = est.fit(unzipped_fit_params.vals);
        std::map<ParamId, double> p_hat_map = zip(p_ids, fr.p_hat);
        std::map<ParamId, double> eta_hat_map = zip(eta_ids, fr.eta_hat);
        std::map<ParamId, double> p_hat_std_map = zip(p_ids, fr.p_hat_std);
        std::map<ParamId, std::map<ParamId, double>> p_hat_corr_map = zip(p_ids, fr.p_hat_correlations);

        FitResultWithMaps out;
        out.p_hat   = std::move(p_hat_map);
        out.eta_hat = std::move(eta_hat_map);
        out.ell_hat = fr.ell_hat;
        out.p_hat_std = std::move(p_hat_std_map);
        out.p_correlations = std::move(p_hat_corr_map);
        out.fit_ok = true;
        this->cache.mle_result = out;

        return out;
    }

    std::set<std::vector<std::pair<double, double>>> confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, CLMethod method = CLMethod::PROJECT);

    void fill_cache() {
        cache.p_specs = this->get_p_specs();
        cache.eta_specs_real = this->get_all_obss_deps();
        for (const auto& [pid, _] : cache.p_specs)
            cache.eta_specs_real.erase(pid);
        
        cache.SigmaEta = this->get_all_correlations();
        cache.exp_obs = this->get_obs_exp();
        cache.SigmaObs = this->get_all_obs_correlations();
    }

    std::map<ParamId, double> get_all_obss_deps() {
        std::unordered_set<ParamId> eta_infos;

        for (const auto& [obsId, qcdOrder] : config.obss) {
            for (auto paramId : obs_int->get_obs_deps(obsId)) {
                if (eta_infos.find(paramId) == eta_infos.end()) {
                    if (pspp->get_param(paramId)->get_combined_std().real() >
                        pspp->get_param(paramId)->get_val() * 1e-6) { // TODO: hardcode à nettoyer
                        eta_infos.insert(paramId);
                    }
                }
            }
        }

        std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);
        std::map<ParamId, double> eta_specs_real_leaf;

        for (auto& paramId : eta_infos_leaf) {
            scalar_t value = pspp->get_param(paramId)->get_val();
            if (pspp->get_param(paramId)->get_combined_std().real() > value.real() * 1e-6) {
                eta_specs_real_leaf[paramId] = value.real();

            } else {
                std::cout << paramId << " does not have real uncertainty" << std::endl;
                std::cout << pspp->get_param(paramId)->get_combined_std().real() << std::endl;
            }
        }
        return eta_specs_real_leaf;
    }

    std::map<ParamId, double> get_p_specs() {
        std::map<ParamId, double> out;
        for (auto elem : config.p_specs) {
            out[elem] = pspp->get_param(elem)->get_val();
        }
        return out;
    }
    std::map<ParamId, std::map<ParamId, double>> get_all_correlations() {
        std::map<ParamId, std::map<ParamId, double>> res;
        CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);
        res = ct.transform(this->cache.eta_specs_real);
        return res;
    }
    
    std::map<ObservableId, std::map<ObservableId, double>> get_all_obs_correlations() {
        std::map<ObservableId, std::map<ObservableId, double>> res;
        CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);

        res = ct.transform(this->cache.exp_obs);
        return res;
    }

    std::map<ObservableId, double> get_obs_exp() {
        std::map<ObservableId, double> out;

        for (auto obs : this->config.obss) {
            out[obs.first] = pspp->get_obs_param(obs.first)->get_val();
        }
        return out;
    }

    void print_cache();

private:
    std::shared_ptr<IModel> obs_int;
    std::shared_ptr<IStatCorrelationProxy> pscp;
    std::shared_ptr<IStatParameterProxy> pspp;
    std::shared_ptr<IStatSourcesProxy> sp;
    StatisticConfig config;
    StatCache cache;
};

#endif