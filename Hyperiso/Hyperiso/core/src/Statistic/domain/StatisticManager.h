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
#include "IStatDependencyPruner.h"
#include "CovarianceTransformer.h"
#include "JointDistribution.h"
#include "RvgNuisanceSampler.h"
#include "MCEngine.h"
#include "DistributionFactory.h"
#include "CopulaFactory.h"
#include "Fit.h"
#include "MarginalConfigFactory.h"

struct StatisticConfig {
    std::vector<ParamId> p_specs;
    std::map<ParamId, MarginalType> override_nuisance_marginals {};
    std::map<BinnedObservableId, MarginalType> override_exp_data_marginals {};
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
    std::map<BinnedObservableId, double> exp_obs;
    std::map<ParamId, double> eta_specs_real;
    std::map<ParamId, double> p_specs;
    std::map<ParamId, std::map<ParamId, double>> SigmaEta;
    std::map<BinnedObservableId, std::map<BinnedObservableId, double>> SigmaObs;

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
        std::shared_ptr<IStatSourcesProxy> sp, std::shared_ptr<IStatDependencyPruner> dp) 
    : config(config), obs_int(obs_int), pscp(pscp), pspp(pspp), sp(sp), dp(dp) {
        obs_int->compute_observables();
    }

    std::vector<std::unique_ptr<IMarginalDistribution>> build_nuisance_marginal_distributions();
    std::unique_ptr<JointDistribution> build_nuisance_distribution();
    std::unique_ptr<JointDistribution> build_exp_data_distribution();

    std::map<BinnedObservableId, GaussianSummary> compute_uncertainties() {
        auto sums = this->compute_uncertainties_and_sampling();

        // Debug print
        for (auto sum : sums.summary) {
            std::cout << sum << std::endl;
        }
        // for (auto _ : unzip(config.obss).ids) {
        //     std::cout << _.str() << std::endl;
        // }

        // for (auto _ : unzip(config.obss).vals) {
        //     std::cout << OrderMapper::str(_) << std::endl;
        // }
        return zip(obs_int->get_obs_ids(), sums.summary);
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
    
    FitResultWithMaps compute_MLE(const std::vector<ParamId>& p_specs) {
        update_cache(std::move(p_specs));
        // Build Likelihood context

        auto unzipped_fit_params = unzip(cache.p_specs);
        auto unzipped_nuisances = unzip(cache.eta_specs_real);
        auto unzipped_exp_obs = unzip(cache.exp_obs);

        std::vector<ParamId> p_ids = unzipped_fit_params.ids;
        std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
        std::vector<BinnedObservableId> obs_ids = unzipped_exp_obs.ids;
        
        LikelihoodContext ctx;
        ctx.nuisance_dist = std::move(build_nuisance_distribution());
        ctx.exp_obs_dist = std::move(build_exp_data_distribution());
        ctx.nuisance_central_values = unzipped_nuisances.vals;
        ctx.exp_obs_values = unzipped_exp_obs.vals;

        // auto model_fn = [this, obs_ids, p_ids, eta_ids] (const Vec& p_vec, const Vec& eta_vec) -> Vec {
        //     auto start = std::chrono::steady_clock::now();
        //     auto pred_map = this->obs_int->predict_optimized(zip(p_ids, p_vec), zip(eta_ids, eta_vec));
        //     auto stop  = std::chrono::steady_clock::now();
        //     auto us = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
        //     // std::cout << "Predict took " << us << " µs" << std::endl;
        //     // for (auto pred : pred_map) {
        //     //     std::cout << ObservableMapper::str(pred.first) << ": ";
        //     //     for (auto ov : pred.second) 
        //     //         std::cout << ov.value << " ";
        //     //     std::cout << std::endl; 
        //     // }
        //     return flatten(pred_map).vals;
        // };

        auto model_fn = [this, obs_ids, p_ids, eta_ids](const Vec& p_vec, const Vec& eta_vec) -> Vec {
            auto pred_map = this->obs_int->predict_optimized(zip(p_ids, p_vec), zip(eta_ids, eta_vec));


            Vec out;
            out.reserve(obs_ids.size());

            for (const auto& bid : obs_ids) {
                // bid.s = ObservableId, bid.p = bin (pair<double,double>)
                const auto& vec = pred_map.at(bid.s);

                // retrouver la bonne entrée dans vec
                // si non binned : bin = {0,0} chez toi, donc match direct
                auto it = std::find_if(vec.begin(), vec.end(), [&](const ObservableValue& ov){
                    auto bin = ov.bin.value_or(std::pair<double,double>{0.,0.});
                    return bin == bid.p; // ou fpeq sur doubles si nécessaire
                });
                if (it == vec.end()) throw std::runtime_error("Missing predicted observable/bin");
                out.push_back(it->value);
            }
            return out;
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

        auto& like = est.like();              // ou stat.get_estimator().like()
        double ell_hat = fr.ell_hat;
        Vector p = fr.p_hat;

        // std::cout << "Scan p2:\n";
        // for (int k=0; k<=40; ++k) {
        //     double p2 = 0.0 + k * 0.01;      // adapte le range !
        //     Vector pp = {p[0], p2};
        //     double d = like.nll_profiled(pp) - ell_hat;
        //     std::cout << p2 << " " << d << "\n";
        // }

        return out;
    }

    std::set<std::vector<std::pair<double, double>>> confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, CLMethod method = CLMethod::PROJECT);

    void update_cache(const std::vector<ParamId>& p_specs = std::vector<ParamId>()) {
        cache.p_specs = this->get_p_specs(p_specs);
        for (const auto& [pid, _] : cache.p_specs)
            dp->detach_parameter(pid.type.value(), pid.block, pid.code);
        cache.eta_specs_real = this->get_all_obss_deps();
        for (const auto& [pid, _] : cache.p_specs)
            cache.eta_specs_real.erase(pid);
        cache.SigmaEta = this->get_all_correlations();
        cache.exp_obs = this->get_obs_exp();
        cache.SigmaObs = this->get_all_obs_correlations();
    }

    std::map<ParamId, double> get_all_obss_deps() {
        std::unordered_set<ParamId> eta_infos;

        for (const auto& obsId : obs_int->get_obs_ids()) {
            for (auto paramId : obs_int->get_obs_deps(obsId.s)) {
                eta_infos.insert(paramId);
            }
        }

        std::unordered_set<ParamId> eta_infos_leaf = this->sp->get_all_leaf_sources(eta_infos);
        std::map<ParamId, double> eta_specs_real_leaf;
        std::map<ParamId, double> delta_rel;

        for (auto& paramId : eta_infos_leaf) {
            double u = pspp->get_param(paramId)->get_combined_std().real();
            if (!std::isfinite(u) || fpeq(std::abs(u), 0.0))
                delta_rel[paramId] = 0;
            else
                delta_rel[paramId] = std::abs(u / pspp->get_param(paramId)->get_val());
        }
            
        
        using T = std::pair<ParamId, double>;
        double delta_rel_max = std::max_element(
            delta_rel.begin(), delta_rel.end(), 
            [] (const T& p, const T& q) { return p.second < q.second; }
        )->second;

        for (auto& paramId : eta_infos_leaf) {
            LOG_INFO("Compared to max relative uncertainty", paramId, delta_rel[paramId] / delta_rel_max);
        }

        double tol = 1e-2; // TODO : make it a parameter
        for (auto& paramId : eta_infos_leaf) {
            if (delta_rel[paramId] / delta_rel_max > tol)
                eta_specs_real_leaf[paramId] = pspp->get_param(paramId)->get_val();
        }

        LOG_INFO("Significant nuisances");
        for (const auto& [pid, val] : eta_specs_real_leaf) {
            LOG_INFO(pid, val);
        }

        return eta_specs_real_leaf;
    }

    std::map<ParamId, double> get_p_specs(const std::vector<ParamId>& p_specs) {
        std::map<ParamId, double> out;
        for (auto elem : p_specs) {            
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
    
    std::map<BinnedObservableId, std::map<BinnedObservableId, double>> get_all_obs_correlations() {
        std::map<BinnedObservableId, std::map<BinnedObservableId, double>> res;
        CovarianceTransformer ct = CovarianceTransformer(pscp, pspp);

        res = ct.transform(this->cache.exp_obs);
        return res;
    }

    std::map<BinnedObservableId, double> get_obs_exp() {
        std::map<BinnedObservableId, double> out;

        for (const auto& obsId : obs_int->get_obs_ids()) {
            out[obsId] = pspp->get_obs_param(obsId)->get_val();
        }
        return out;
    }

    void print_cache();

private:
    std::shared_ptr<IModel> obs_int;
    std::shared_ptr<IStatCorrelationProxy> pscp;
    std::shared_ptr<IStatParameterProxy> pspp;
    std::shared_ptr<IStatSourcesProxy> sp;
    std::shared_ptr<IStatDependencyPruner> dp;
    StatisticConfig config;
    StatCache cache;
};

#endif