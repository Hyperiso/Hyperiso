#include <vector>
#include "Include.h"
#include "IModel.h"
#include "IStatCorrelationProxy.h"
#include "IStatParameterProxy.h"
#include "CovarianceTransformer.h"
#include "RandomVectorGenerator.h"
#include "RvgNuisanceSampler.h"
#include "MonteCarloPredictorGeneric.h"
#include "DistributionFactory.h"

struct StatisticConfig {
    std::map<ObservableId, QCDOrder> obss;
    std::vector<ParamId> p_specs;
};

struct StatCache {
    std::map<ObservableId, double> exp_obs;
    std::map<ParamId, double> eta_specs_real;
    std::map<ParamId, double> p_specs;
    std::map<ParamId, std::map<ParamId, double>> SigmaEta;
    std::map<ObservableId, std::map<ObservableId, double>> SigmaObs;
};

class StatisticManager {
public:
    StatisticManager(StatisticConfig config, std::shared_ptr<IModel> obs_int, std::shared_ptr<IStatCorrelationProxy> pscp, std::shared_ptr<IStatParameterProxy> pspp) : config(config) {
        obs_int->add_observables(config.obss);

    }

    std::map<ObservableId, double> compute_uncertainties() {
        unsigned int seed = std::random_device{}();
        auto dist = DistributionFactory::create("gaussian", seed);
        auto decomp = std::make_unique<CholeskyDecomposition>();

        RandomVectorGenerator rvg(std::move(dist), std::move(decomp));

        std::cout << "RandomVectorGenerator created" << std::endl;

        RvgNuisanceSampler sampler(rvg);

        std::cout << "Creating MonteCarloPredictor" << std::endl;

        // MC prediction with pluggable sampler
        MonteCarloPredictor2 mc(this->obs_int, sampler, cache.eta_specs_real, cache.SigmaEta, {10, 0.2});
        std::mt19937 rng(1234);

        std::cout << "MonteCarloPredictor created" << std::endl;

        // Vec p_test{-4.5, 0.0};
        auto sums = mc.summarize(this->cache.p_specs, rng);
    }
    void fill_cache() {
        cache.eta_specs_real = this->get_all_obss_deps();
        cache.SigmaEta = this->get_all_correlations();
        cache.exp_obs = this->get_obs_exp();
        cache.SigmaObs = this->get_all_obs_correlations();
        cache.p_specs = this->get_p_specs();
    }

    std::map<ParamId, double> get_all_obss_deps() {
        std::map<ParamId, double> eta_specs_real;
        Vec eta_mean_real;
        for (auto elem : config.obss) {
            for (auto _ : obs_int->get_obs_deps(elem.first))
            if (!(std::find(eta_specs_real.begin(), eta_specs_real.end(), _) != eta_specs_real.end())) {
                if (pspp->get_param(_)->get_combined_std().real() > pspp->get_param(_)->get_val() *1e-6) { //TODO bad harcoded
                    eta_specs_real[_] = pspp->get_param(_)->get_val();
                }
            }
        }
        return eta_specs_real;
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
    StatisticConfig config;
    StatCache cache;
};