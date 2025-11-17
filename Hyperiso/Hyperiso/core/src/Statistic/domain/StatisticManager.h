#include <vector>
#include "Include.h"
#include "IModel.h"
#include "IStatCorrelationProxy.h"
#include "IStatParameterProxy.h"

struct StatisticConfig {
    std::map<ObservableId, QCDOrder> obss;
};

struct StatCache {
    std::map<ObservableId, double> exp_obs;
    std::map<ParamId, double> eta_specs_real;
};

class StatisticManager {
public:
    StatisticManager(StatisticConfig config, std::shared_ptr<IModel> obs_int, std::shared_ptr<IStatCorrelationProxy> pscp, std::shared_ptr<IStatParameterProxy> pspp) : config(config) {
        obs_int->add_observables(config.obss);

    }

    std::map<ParamId, double> get_all_obss_deps() {
        std::map<ParamId, double> eta_specs_real;
        Vec eta_mean_real;
        for (auto elem : config.obss) {
            for (auto _ : obs_int->get_obs_deps(elem.first))
            if (!(std::find(eta_specs_real.begin(), eta_specs_real.end(), _) != eta_specs_real.end())) {
                eta_specs_real[_] = pspp->get_param(_)->get_val();
            }
        }
        return eta_specs_real;
    }
private:
    std::shared_ptr<IModel> obs_int;
    std::shared_ptr<IStatCorrelationProxy> pscp;
    std::shared_ptr<IStatParameterProxy> pspp;
    StatisticConfig config;
    StatCache cache;
};