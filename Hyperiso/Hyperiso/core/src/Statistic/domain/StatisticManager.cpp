#include "StatisticManager.h"

std::unique_ptr<JointDistribution> StatisticManager::build_nuisance_distribution() {
    // unsigned int seed = std::random_device{}();
    unsigned int seed = 123456u; //TODO : Niels ta faute
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

    std::unique_ptr<ICopula> copula;
    if (config.nuisance_copula_type == CopulaType::GAUSSIAN) {
        GaussianCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    } else if (config.nuisance_copula_type == CopulaType::STUDENT_T) {
        StudentTCopulaConfig copula_cfg;
        copula_cfg.R = RealMatrix(unzip(cache.SigmaEta).vals);
        copula_cfg.nu = cache.eta_specs_real.size() - 1;
        copula = CopulaFactory::create(config.nuisance_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));
}

std::unique_ptr<JointDistribution> StatisticManager::build_exp_data_distribution() {
    unsigned int seed = std::random_device{}();
    std::map<ObservableId, MarginalType> exp_data_marginals;

    for (auto& [oid, v] : cache.exp_obs) {
        exp_data_marginals.emplace(oid, MarginalType::GAUSSIAN);
    }

    for (auto& [oid, mt] : config.override_exp_data_marginals) {
        if (!exp_data_marginals.contains(oid)) {
            LOG_WARN("Observable", ObservableMapper::str(oid), "is not a taken into account in this fit.");
            continue;
        }
        exp_data_marginals.at(oid) = mt;
    }

    auto unzipped = unzip(exp_data_marginals);
    std::vector<ObservableId> obs_ids = unzipped.ids;
    std::vector<std::unique_ptr<IMarginalDistribution>> marginals;

    for (auto& [oid, mt] : exp_data_marginals) {
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
        copula_cfg.nu = config.obss.size() - 1;
        copula = CopulaFactory::create(config.exp_data_copula_type, copula_cfg, seed);
    }

    return std::make_unique<JointDistribution>(std::move(marginals), std::move(copula));
}
