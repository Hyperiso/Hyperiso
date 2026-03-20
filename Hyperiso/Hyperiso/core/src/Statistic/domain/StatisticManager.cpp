#include "StatisticManager.h"

std::vector<std::unique_ptr<IMarginalDistribution>> StatisticManager::build_nuisance_marginal_distributions() {
    // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
    unsigned int seed = 123456u;
    std::map<ParamId, MarginalType> nuisance_marginals;

    for (auto& [pid, v] : cache.eta_specs_real) {
        nuisance_marginals.emplace(pid, MarginalType::GAUSSIAN);
    }

    for (auto& [pid, mt] : config.override_nuisance_marginals) {
        for (auto elem : nuisance_marginals) {
            std::cout << elem.first << std::endl;
            std::cout << ParameterTypeMapper::str(elem.first.type.value_or(ParameterType::WILSON)) << std::endl;
        }
        if (!nuisance_marginals.contains(pid)) {
            LOG_WARN("Parameter", pid, "is not a nuisance for the selected observables");
            continue;
        }
        nuisance_marginals.at(pid) = mt;
    }

    auto unzipped = unzip(nuisance_marginals);
    std::vector<ParamId> nuisance_ids = unzipped.ids;
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

    for (auto& [pid, mt] : nuisance_marginals) {
        MarginalConfig cfg = MarginalConfigFactory().create(pid, mt);
        auto m_ptr = DistributionFactory::create(mt, cfg, seed);
        marginals.emplace_back(std::move(m_ptr));
    }

    return marginals;
}

std::unique_ptr<JointDistribution> StatisticManager::build_nuisance_distribution() {
    // unsigned int seed = std::random_device{}(); TODO : investigate how random seed leads to unstable results
    unsigned int seed = 123456u;

    std::unique_ptr<ICopula> copula;
    if (config.nuisance_copula_type == CopulaType::GAUSSIAN) {
        GaussianCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    } else if (config.nuisance_copula_type == CopulaType::STUDENT_T) {
        StudentTCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula_cfg.nu = cache.eta_specs_real.size();
        copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(
        std::move(build_nuisance_marginal_distributions()), 
        std::move(copula)
    );
}

std::unique_ptr<JointDistribution> StatisticManager::build_exp_data_distribution() {
    unsigned int seed = std::random_device{}();
    std::map<ExperimentObs, MarginalType> exp_data_marginals;

    for (auto& [oid, v] : cache.exp_obs) {
        exp_data_marginals.emplace(oid, MarginalType::GAUSSIAN);
    }

    for (auto& [oid, mt] : config.override_exp_data_marginals) {
        if (!exp_data_marginals.contains(oid)) {
            LOG_WARN("Observable", ObservableMapper::str(oid.obs.s), "is not a taken into account in this fit.");
            continue;
        }
        exp_data_marginals.at(oid) = mt;
    }

    auto unzipped = unzip(exp_data_marginals);
    std::vector<ExperimentObs> obs_ids = unzipped.ids;
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

    for (auto& [oid, mt] : exp_data_marginals) {//TODO : checkkkkkk
        MarginalConfig cfg = MarginalConfigFactory().create(oid, mt);
        auto m_ptr = DistributionFactory::create(mt, cfg, seed);
        marginals.emplace_back(std::move(m_ptr));
    }

    std::unique_ptr<ICopula> copula;
    if (config.exp_data_copula_type == CopulaType::GAUSSIAN) {
        GaussianCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
        copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
    } else if (config.exp_data_copula_type == CopulaType::STUDENT_T) {
        StudentTCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaObs).vals);
        copula_cfg.nu = obs_int->n_observables() - 1;
        copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));
}

std::set<std::vector<std::pair<double, double>>> StatisticManager::confidence_contour(ParamId p1, ParamId p2, double z, std::array<double, 4> bounds, CLMethod method) {
    // if (!(this->cache.p_specs.contains(p1) && this->cache.p_specs.contains(p2)))
    //     throw std::invalid_argument("Invalid parameter for confidence level.");

    // if (!this->cache.mle_result.fit_ok)
    //     throw std::runtime_error("Please perform MLE fit before extracting contours");

    // // Build joint distribution (gaussian marginals + gaussian copula) for the remaining fit parameters
    // unsigned int seed = 123456u;

    // // Marginals
    // std::vector<std::unique_ptr<IMarginalDistribution>> marginals;
    // for (auto& [p, _] : cache.p_specs) {
    //     if (p == p1 || p == p2) continue;

    //     GaussianMarginalCfg cfg(cache.mle_result.p_hat.at(p), cache.mle_result.p_hat_std.at(p));
    //     auto m_ptr = DistributionFactory::create(MarginalType::GAUSSIAN, cfg, seed);
    //     marginals.emplace_back(std::move(m_ptr));
    // }

    // auto nuisance_marginals = build_nuisance_marginal_distributions();
    // marginals.insert(
    //     marginals.end(),
    //     std::make_move_iterator(nuisance_marginals.begin()),
    //     std::make_move_iterator(nuisance_marginals.end())
    // );

    // // Copula
    // auto p_ids = unzip(cache.p_specs).ids;
    // RealMatrix R = unzip(cache.mle_result.p_correlations).vals;
    // std::vector<std::vector<ParamId>::iterator> to_erase;
    // for (auto& p : {p1, p2}) {
    //     std::ptrdiff_t idx = std::distance(p_ids.begin(), std::find(p_ids.begin(), p_ids.end(), p));
    //     to_erase.emplace_back(p_ids.begin() + idx);
    //     R.remove_row_and_column(idx);
    // }

    // for (auto& it : to_erase)
    //     p_ids.erase(it);

    // RealMatrix combined_R = block_diag(RealMatrix(unzip(cache.SigmaEta).vals), R);

    // std::unique_ptr<ICopula> copula;
    // if (config.nuisance_copula_type == CopulaType::GAUSSIAN) {
    //     GaussianCopulaConfig copula_cfg;
    //     copula_cfg.R = combined_R;
    //     copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    // } else if (config.nuisance_copula_type == CopulaType::STUDENT_T) {
    //     StudentTCopulaConfig copula_cfg;
    //     copula_cfg.R = combined_R;
    //     copula_cfg.nu = combined_R.rows();
    //     copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    // }

    // auto combined_distribution = std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));

    // Dispatch ids
    // std::vector<ParamId> p_ids_2 {p1, p2};
    // auto unzipped_nuisances = unzip(cache.eta_specs_real);
    // auto unzipped_exp_obs = unzip(cache.exp_obs);
    // std::vector<ParamId> eta_ids = unzipped_nuisances.ids;
    // std::vector<ExperimentObs> obs_ids = unzipped_exp_obs.ids; 
    // eta_ids.insert(eta_ids.end(), p_ids.begin(), p_ids.end());

    // auto model_fn = [this, obs_ids, p_ids, eta_ids] (const Vec& p_vec, const Vec& eta_vec) -> Vec {
    //     auto pred_map = this->obs_int->predict_optimized(zip(p_ids, p_vec), zip(eta_ids, eta_vec));
    //     return flatten(pred_map).vals;
    // };

    // // Create likelihood
    // LikelihoodContext ctx;
    // ctx.exp_obs_dist = std::move(build_exp_data_distribution());
    // ctx.nuisance_dist = std::move(combined_distribution);
    // ctx.exp_obs_values = unzipped_exp_obs.vals;
    // ctx.nuisance_central_values = unzipped_nuisances.vals;

    // for (auto& pid : p_ids)
    //     ctx.nuisance_central_values.emplace_back(cache.mle_result.p_hat.at(pid));

    // MLEstimator fitter(std::move(ctx), model_fn);
    // return fitter.contour(z, bounds, this->cache.mle_result.ell_hat);
}

void StatisticManager::print_cache()
{
    for (auto elem : cache.eta_specs_real) {
        std::cout << " eta_specs_real : " << elem.first << " = " << elem.second;
    }
    std::cout << std::endl;

    for (auto elem : cache.SigmaEta) {
        for (auto elem2 : elem.second) {
            std::cout << " SigmaEta : " << elem.first << " | " << elem2.first << " = " << elem2.second;
        }
        std::cout << std::endl;
    }

    for (auto elem : cache.exp_obs) {
        std::cout << " exp_obs : " << elem.first.str() << " = " << elem.second;
    }
    std::cout << std::endl;

    for (auto elem : cache.SigmaObs) {
        for (auto elem2 : elem.second) {
            std::cout << " SigmaObs : " << elem.first.str() << " | " << elem2.first.str() << " = " << elem2.second;
        }
        std::cout << std::endl;
    }

    for (auto elem : cache.p_specs) {
        std::cout << " p_specs : " << elem.first << " = " << elem.second;
    }
    std::cout << std::endl;

    std::cout << "eta size : " << this->cache.eta_specs_real.size() << std::endl;
    std::cout << "etasigma size : " << this->cache.eta_specs_real.size() << " | " << this->cache.SigmaEta.begin()->second.size() << std::endl;
    std::cout << "p_specs size : " << this->cache.p_specs.size() << std::endl;
    std::cout << "exp_obs : " << this->cache.exp_obs.size() << std::endl;
    std::cout << "SigmaObs size : " << this->cache.SigmaObs.size() << " | " << this->cache.SigmaObs.begin()->second.size() << std::endl;
}
