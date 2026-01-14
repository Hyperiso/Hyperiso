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

    std::size_t MLE_max_iter = 100;
    double MLE_tol = 1e-8;
};

struct FitResultWithMaps {
    std::map<ParamId, double> p_hat;
    std::map<ParamId, double> eta_hat;
    double ell_hat{0.0};
};

struct StatCache {
    std::map<ObservableId, double> exp_obs;
    std::map<ParamId, double> eta_specs_real;
    std::map<ParamId, double> p_specs;
    std::map<ParamId, std::map<ParamId, double>> SigmaEta;
    std::map<ObservableId, std::map<ObservableId, double>> SigmaObs;

    FitResultWithMaps mle_result;
};

class StatisticManager {
public:
    StatisticManager(StatisticConfig config, std::shared_ptr<IModel> obs_int, 
        std::shared_ptr<IStatCorrelationProxy> pscp, std::shared_ptr<IStatParameterProxy> pspp,
        std::shared_ptr<IStatSourcesProxy> sp) 
    : config(config), obs_int(obs_int), pscp(pscp), pspp(pspp), sp(sp) {
        obs_int->add_observables(config.obss);
    }

    std::unique_ptr<JointDistribution> build_nuisance_distribution();
    std::unique_ptr<JointDistribution> build_exp_data_distribution();

    // std::pair<double,double> compute_CI_1d_95(
    //     const ParamId& pid_to_scan,
    //     double p_min, double p_max,
    //     int grid_points
    // ) {

    //     std::vector<ObservableId> obs_order;
    //     obs_order.reserve(cache.exp_obs.size());
    //     for (const auto& [oid, val] : cache.exp_obs) {
    //         obs_order.push_back(oid);
    //     }

    //     std::vector<ParamId> eta_order;
    //     eta_order.reserve(cache.eta_specs_real.size());
    //     for (const auto& [pid, val] : cache.eta_specs_real) {
    //         eta_order.push_back(pid);
    //     }

    //     std::vector<ParamId> p_order;
    //     p_order.reserve(cache.p_specs.size());
    //     for (const auto& [pid, val] : cache.p_specs) {
    //         p_order.push_back(pid);
    //     }

    //     auto vec_from_param_map = [](const std::map<ParamId, double>& m,
    //                                 const std::vector<ParamId>& order) {
    //         Vec v(order.size());
    //         for (std::size_t i = 0; i < order.size(); ++i) {
    //             auto it = m.find(order[i]);
    //             if (it == m.end())
    //                 throw std::runtime_error("vec_from_param_map: missing ParamId in map");
    //             v[i] = it->second;
    //         }
    //         return v;
    //     };

    //     auto vec_from_obs_map = [](const std::map<ObservableId, double>& m,
    //                             const std::vector<ObservableId>& order) {
    //         Vec v(order.size());
    //         for (std::size_t i = 0; i < order.size(); ++i) {
    //             auto it = m.find(order[i]);
    //             if (it == m.end())
    //                 throw std::runtime_error("vec_from_obs_map: missing ObservableId in map");
    //             v[i] = it->second;
    //         }
    //         return v;
    //     };

    //     auto mat_from_eta_cov_map =
    //         [&](const std::map<ParamId, std::map<ParamId, double>>& cov,
    //             const std::vector<ParamId>& order) {
    //             Matrix M(order.size(), Vec(order.size(), 0.0));
    //             for (std::size_t i = 0; i < order.size(); ++i) {
    //                 for (std::size_t j = 0; j < order.size(); ++j) {
    //                     auto it_row = cov.find(order[i]);
    //                     if (it_row != cov.end()) {
    //                         auto it_col = it_row->second.find(order[j]);
    //                         if (it_col != it_row->second.end()) {
    //                             M[i][j] = it_col->second;
    //                         } else {
    //                             M[i][j] = 0.0;
    //                         }
    //                     } else {
    //                         M[i][j] = 0.0;
    //                     }
    //                 }
    //             }
    //             return M;
    //         };

    //     auto mat_from_obs_cov_map =
    //         [&](const std::map<ObservableId, std::map<ObservableId, double>>& cov,
    //             const std::vector<ObservableId>& order) {
    //             Matrix M(order.size(), Vec(order.size(), 0.0));
    //             for (std::size_t i = 0; i < order.size(); ++i) {
    //                 for (std::size_t j = 0; j < order.size(); ++j) {
    //                     auto it_row = cov.find(order[i]);
    //                     if (it_row != cov.end()) {
    //                         auto it_col = it_row->second.find(order[j]);
    //                         if (it_col != it_row->second.end()) {
    //                             M[i][j] = it_col->second;
    //                         } else {
    //                             M[i][j] = 0.0;
    //                         }
    //                     } else {
    //                         M[i][j] = 0.0;
    //                     }
    //                 }
    //             }
    //             return M;
    //         };


    //     Vec Oexp_vec    = vec_from_obs_map(cache.exp_obs,          obs_order);
    //     Vec eta_bar_vec = vec_from_param_map(cache.eta_specs_real, eta_order);

    //     Matrix SigmaO_mat   = mat_from_obs_cov_map(cache.SigmaObs,  obs_order);
    //     Matrix SigmaEta_mat = mat_from_eta_cov_map(cache.SigmaEta,  eta_order);

    //     SPDMatrix SO = SPDMatrix::cholesky(SigmaO_mat);
    //     SPDMatrix SE = SPDMatrix::cholesky(SigmaEta_mat);

    //     LikelihoodContext ctx{Oexp_vec, SO, eta_bar_vec, SE};


    //     auto model_fn = [this, obs_order, p_order, eta_order]
    //                     (const Vec& p_vec, const Vec& eta_vec) -> Vec {
    //         std::map<ParamId, double> p_map;
    //         for (std::size_t i = 0; i < p_order.size(); ++i) {
    //             p_map[p_order[i]] = p_vec[i];
    //         }

    //         std::map<ParamId, double> eta_map;
    //         for (std::size_t i = 0; i < eta_order.size(); ++i) {
    //             eta_map[eta_order[i]] = eta_vec[i];
    //         }

    //         auto pred_map = this->obs_int->predict_optimized(p_map, eta_map);

    //         Vec pred_vec(obs_order.size());
    //         for (std::size_t i = 0; i < obs_order.size(); ++i) {
    //             auto it = pred_map.find(obs_order[i]);
    //             if (it == pred_map.end()) {
    //                 throw std::runtime_error("model_fn: missing prediction for observable");
    //             }
    //             pred_vec[i] = it->second;
    //         }

    //         return pred_vec;
    //     };

    //     MLEstimator est(ctx, model_fn, this->config.MLE_max_iter, this->config.MLE_tol);


    //     const FitResultWithMaps& fr_map = this->cache.mle_result;

    //     FitResult fr_vec;
    //     fr_vec.p_hat   = vec_from_param_map(fr_map.p_hat,   p_order);
    //     fr_vec.eta_hat = vec_from_param_map(fr_map.eta_hat, eta_order);
    //     fr_vec.ell_hat = fr_map.ell_hat;

    //     std::cout << "MLE fit done (from cache), ell_hat = " << fr_vec.ell_hat << std::endl;

    //     std::cout << "finding thr95 " << std::endl;
    //     const double thr95 = gsl_cdf_chisq_Pinv(0.95, 1);
    //     std::cout << "thr95 : " << thr95 << std::endl;


    //     int idx_scan = -1;
    //     for (std::size_t i = 0; i < p_order.size(); ++i) {
    //         if (p_order[i] == pid_to_scan) {
    //             idx_scan = static_cast<int>(i);
    //             break;
    //         }
    //     }
    //     if (idx_scan < 0) {
    //         throw std::runtime_error("compute_CI_1d_95: pid_to_scan not found in p_order");
    //     }


    //     Vec eta_init = eta_bar_vec;

    //     auto T = [&](double p_scan_value) {
    //         Vec p_test = fr_vec.p_hat;  
    //         p_test[idx_scan] = p_scan_value; 
    //         return est.wilks_T(p_test, fr_vec, eta_init);
    //     };


    //     double a = p_min;
    //     double b = p_max;
    //     int N = grid_points;

    //     auto find_crossing = [&](double lo, double hi, int N, bool want_left) -> double {

    //     double prev = want_left ? hi : lo;
    //     double prevT = T(prev);

    //     for (int i = 1; i <= N; ++i) {
    //         double x;
    //         if (want_left) {
    //             x = hi - (hi - lo) * i / double(N);
    //         } else {
    //             x = lo + (hi - lo) * i / double(N);
    //         }

    //         double t = T(x);

    //         if ((prevT - thr95) * (t - thr95) <= 0.0) {
    //             return x;
    //         }

    //         prev = x;
    //         prevT = t;
    //     }

    //     return std::nan("");
    // };

    //     // double left  = std::nan("");
    //     // double right = std::nan("");

    //     // double prev  = a;
    //     // double prevT = T(prev);

    //     // for (int i = 1; i <= N; ++i) {
    //     //     double x = a + (b - a) * i / double(N);
    //     //     double t = T(x);

    //     //     if (std::isnan(left) && (prevT - thr95) * (t - thr95) <= 0.0)
    //     //         left = prev;

    //     //     if ((prevT - thr95) * (t - thr95) <= 0.0)
    //     //         right = x;

    //     //     prev  = x;
    //     //     prevT = t;
    //     //     std::cout << "scan p=" << x << "  T(p)=" << t << "\n";
    //     // }

    //     double phat = fr_vec.p_hat[idx_scan];
        
    //     double left  = find_crossing(p_min, phat, N, true);
    //     double right = find_crossing(phat, p_max, N, false);

    //     std::cout << "95% CI for param " << pid_to_scan
    //             << ": [" << left << "," << right << "]" << std::endl;

    //     return {left, right};
    // }

    std::map<ObservableId, GaussianSummary> compute_uncertainties() {
        auto rvg = build_nuisance_distribution();
        std::vector<ParamId> nuisance_ids = unzip(cache.eta_specs_real).ids;
        RvgNuisanceSampler sampler(nuisance_ids, std::move(rvg));
        MonteCarloEngine mc(this->obs_int, sampler, {this->config.MC_draws, this->config.skew_abs_threshold});

        auto sums = mc.summarize(this->cache.p_specs);

        // Debug print
        for (auto sum : sums) {
            std::cout << sum << std::endl;
        }

        return zip(unzip(config.obss).ids, sums);
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

            Vec pred_vec(obs_ids.size());
            for (std::size_t i = 0; i < obs_ids.size(); ++i) {
                if (!pred_map.contains(obs_ids[i]))
                    throw std::runtime_error("model_fn: missing prediction for observable");
                
                pred_vec[i] = pred_map[obs_ids[i]];
            }

            return pred_vec;
        };

        MLEstimator est(std::move(ctx), model_fn, this->config.MLE_max_iter, this->config.MLE_tol);

        std::cout << "Now doing MLE : " << std::endl;

        /*
        auto debug_scan_around_start = [&](MLEstimator& est,
                                    const Vec& p0, const Vec& eta0,
                                    const std::vector<ParamId>& p_order,
                                    const std::vector<ParamId>& eta_order) {
            std::cout << "\n=== DEBUG scan around starting point ===\n";
            double ell0 = est.like().ell(p0, eta0);
            std::cout << "ell(p0, eta0) = " << ell0 << "\n";
            std::cout << std::setprecision(17);
            // scan sur les p
            for (std::size_t i = 0; i < p0.size(); ++i) {
                Vec p_plus = p0;
                Vec p_minus = p0;
                double delta = 0.1 * std::abs(p0[i]); // par ex 10% de la valeur
                // if (delta == 0.0) delta = 1e-3;       // fallback
                double floor = 1000.0 * std::numeric_limits<double>::epsilon() * (std::abs(p0[i]) + 1.0);
                delta = std::max(delta, floor);
                p_plus[i]  += delta;
                p_minus[i] -= delta;

                double ell_plus  = est.like().ell(p_plus, eta0);
                double ell_minus = est.like().ell(p_minus, eta0);

                std::cout << "Param " << p_order[i]
                        << " : ell(+delta)=" << ell_plus
                        << ", ell(-delta)=" << ell_minus << "\n";
            }

            // scan sur les eta
            for (std::size_t j = 0; j < eta0.size(); ++j) {
                Vec eta_plus = eta0;
                Vec eta_minus = eta0;
                double delta = 0.1 * std::abs(eta0[j]);
                if (delta == 0.0) delta = 1e-3;

                eta_plus[j]  += delta;
                eta_minus[j] -= delta;

                double ell_plus  = est.like().ell(p0, eta_plus);
                double ell_minus = est.like().ell(p0, eta_minus);

                std::cout << "Eta " << eta_order[j]
                        << " : ell(+delta)=" << ell_plus
                        << ", ell(-delta)=" << ell_minus << "\n";
            }
        };
        debug_scan_around_start(est, p0, eta0, p_order, eta_order);
        */

        FitResult fr = est.fit(unzipped_fit_params.vals);

        std::cout << "MLE fit done: ell_hat = " << std::setprecision(5) << fr.ell_hat << std::endl;

        std::map<ParamId, double> p_hat_map = zip(p_ids, fr.p_hat);
        std::map<ParamId, double> eta_hat_map = zip(eta_ids, fr.eta_hat);

        for (const auto& [pid, val] : p_hat_map) {
            std::cout << "p_hat[" << pid << "] = " << val << std::endl;
        }
        for (const auto& [pid, val] : eta_hat_map) {
            std::cout << "eta_hat[" << pid << "] = " << val << std::endl;
        }

        FitResultWithMaps out;
        out.p_hat   = std::move(p_hat_map);
        out.eta_hat = std::move(eta_hat_map);
        out.ell_hat = fr.ell_hat;

        this->cache.mle_result = out;

        return out;
    }

    void fill_cache() {

        cache.p_specs = this->get_p_specs();
        cache.eta_specs_real = this->get_all_obss_deps();
        for (const auto& [pid, _] : cache.p_specs) {
            cache.eta_specs_real.erase(pid);
        }
        for (auto elem : cache.eta_specs_real) {
            std::cout << " eta_specs_real : " << elem.first << " = " << elem.second;
        }
        std::cout << std::endl;
        cache.SigmaEta = this->get_all_correlations();

        for (auto elem : cache.SigmaEta) {
            for (auto elem2 : elem.second) {
                std::cout << " SigmaEta : " << elem.first << " | " << elem2.first << " = " << elem2.second;
            }
            std::cout << std::endl;
        }
        cache.exp_obs = this->get_obs_exp();
        for (auto elem : cache.exp_obs) {
            std::cout << " exp_obs : " << elem.first.str() << " = " << elem.second;
        }

        std::cout << std::endl;
        cache.SigmaObs = this->get_all_obs_correlations();

        for (auto elem : cache.SigmaObs) {
            for (auto elem2 : elem.second) {
                std::cout << " SigmaObs : " << elem.first.str() << " | " << elem2.first.str() << " = " << elem2.second;
            }
            std::cout << std::endl;
        }

        // cache.p_specs = this->get_p_specs();

        for (auto elem : cache.p_specs) {
            std::cout << " p_specs : " << elem.first << " = " << elem.second;
        }

        std::cout << std::endl;
        std::cout << "eta size : " << this->cache.eta_specs_real.size() << std::endl;
        std::cout << "etasigma size : " << this->cache.eta_specs_real.size() << " | " << this->cache.SigmaEta.at(ParamId(ParameterType::SM, "VCKMIN", 1)).size() << std::endl;
        std::cout << "p_specs size : " << this->cache.p_specs.size() << std::endl;
        std::cout << "exp_obs : " << this->cache.exp_obs.size() << std::endl;
        std::cout << "SigmaObs size : " << this->cache.SigmaObs.size() << " | " << this->cache.SigmaObs[ObservableMapper::to_id(Observables::BR_BD_MUMU)].size() << std::endl;
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

private:
    std::shared_ptr<IModel> obs_int;
    std::shared_ptr<IStatCorrelationProxy> pscp;
    std::shared_ptr<IStatParameterProxy> pspp;
    std::shared_ptr<IStatSourcesProxy> sp;
    StatisticConfig config;
    StatCache cache;
};

#endif