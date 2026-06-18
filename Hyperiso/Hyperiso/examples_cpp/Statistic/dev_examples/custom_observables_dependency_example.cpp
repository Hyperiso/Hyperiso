#include <complex>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ObservableInterface.h"
#include "ObservableInterfaceProxy.h"
#include "StatParamOptimizerProxy.h"
#include "StatisticInterface.h"
#include "mapper_hub.hpp"

namespace {

void print_param_set(const std::string& title, const std::unordered_set<ParamId>& deps) {
    std::cout << "\n" << title << "\n";
    if (deps.empty()) {
        std::cout << "  <empty>\n";
        return;
    }
    for (const auto& pid : deps) {
        std::cout << "  " << pid << "\n";
    }
}

void print_param_map(const std::string& title, const std::map<ParamId, double>& deps) {
    std::cout << "\n" << title << "\n";
    if (deps.empty()) {
        std::cout << "  <empty>\n";
        return;
    }
    for (const auto& [pid, value] : deps) {
        std::cout << "  " << pid << " = " << value << "\n";
    }
}

void print_predictions(const std::string& title,
                       const std::map<ObservableId, std::vector<ObservableValue>>& predictions) {
    std::cout << "\n" << title << "\n";
    for (const auto& [obs, values] : predictions) {
        for (const auto& value : values) {
            std::cout << "  " << ObservableMapper::str(obs) << " = " << value.value;
            if (value.bin.has_value()) {
                std::cout << " in [" << value.bin->first << ", " << value.bin->second << "]";
            }
            std::cout << "\n";
        }
    }
}

} // namespace

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    HyperisoConfig config_hyp;
    config_hyp.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config_hyp);
    init_all_builtins();

    auto oi = std::make_shared<ObservableInterface>();

    // Three ordinary input parameters used by the custom Wilson lambdas and observables.
    // These are deliberately standard Hyperiso ParamId values so the Statistic layer can
    // find uncertainties/correlations and can mutate them through StatParamOptimizerProxy.
    const ParamId f_bd(ParameterType::FLAVOR, "FCONST", LhaID(511, 1));
    const ParamId f_bs(ParameterType::FLAVOR, "FCONST", LhaID(531, 1));
    const ParamId m_top(ParameterType::SM, "SMINPUTS", LhaID(6));

    // Custom Wilson group. The dependency sources declared in set_matching(...) are
    // automatically propagated to every observable in the LambdaDecayConfig below.
    GroupMapper::register_custom("STAT_DEV_WILSON", {"stat-dev-wilson"});
    WCoefMapper::register_custom("C_STAT_A", {"cstatA"}, std::pair<int, int>{993001, 1});
    WCoefMapper::register_custom("C_STAT_B", {"cstatB"}, std::pair<int, int>{993001, 2});

    WGroupId group = GroupMapper::id_of("stat-dev-wilson");
    WCoefId c_a = WCoefMapper::id_of("cstatA");
    WCoefId c_b = WCoefMapper::id_of("cstatB");

    CustomWilsonCoefficientConfig wc_a(c_a);
    wc_a.set_matching(
        QCDOrder::LO,
        {f_bd, m_top},
        [f_bd, m_top](const ParamSrc& src) {
            const double f = std::real(src.get_val(f_bd));
            const double mt = std::real(src.get_val(m_top));
            return 0.2 + 0.5 * f + 1.0e-3 * (mt - 172.0);
        },
        ContributionType::SM
    );

    CustomWilsonCoefficientConfig wc_b(c_b);
    wc_b.set_matching(
        QCDOrder::LO,
        {f_bs, m_top},
        [f_bs, m_top](const ParamSrc& src) {
            const double f = std::real(src.get_val(f_bs));
            const double mt = std::real(src.get_val(m_top));
            return -0.1 + 0.3 * f - 5.0e-4 * (mt - 172.0);
        },
        ContributionType::SM
    );

    CustomWilsonGroupConfig wc_group(group);
    wc_group.matching_scale = 81.0;
    wc_group.hadronic_scale = 4.8;
    wc_group.order = QCDOrder::LO;
    wc_group.contribution = ContributionType::SM;
    wc_group.add_coefficient(wc_a).add_coefficient(wc_b);

    // Optional running lambda. Identity running would also work, but this makes the
    // example explicit and shows that WCoefId custom coefficients flow through running.
    wc_group.set_running(
        WilsonBasis::B_STANDARD,
        QCDOrder::LO,
        {},
        [c_a, c_b](const std::unordered_map<QCDOrder, std::unordered_map<WCoefId, scalar_t>>& matching,
                   const BlockSrc&) {
            std::unordered_map<WCoefId, scalar_t> out;
            out[c_a] = 0.95 * matching.at(QCDOrder::LO).at(c_a);
            out[c_b] = 1.05 * matching.at(QCDOrder::LO).at(c_b);
            return out;
        }
    );

    // Custom decay with two observables. Each observable can add its own direct
    // dependencies. The custom-Wilson matching dependencies above are additionally
    // propagated because LambdaDecayConfig::propagate_custom_wilson_dependencies=true.
    LambdaDecayConfig decay;
    decay.canonical = "STAT_DEV_DECAY";
    decay.aliases = {"stat-dev-decay"};
    decay.matching_scale = 81.0;
    decay.hadronic_scale = 4.8;
    decay.order = QCDOrder::LO;
    decay.custom_wilson_groups.push_back(wc_group);

    LambdaObservableConfig obs_a = LambdaObservableConfig::scalar(
        "STAT_DEV_OBS_A",
        [group, c_a, f_bd](LambdaDecay& ctx, ObservableId id) {
            (void)id;
            const double f = std::real(ctx.FLAVOR()(f_bd, DataType::VALUE));
            const double c = std::real(ctx.W().getFR(group, c_a, QCDOrder::LO, ContributionType::SM));
            return c + f;
        }
    );
    obs_a.aliases = {"stat-dev-obs-a"};
    obs_a.flha = LhaID(994001, 1);
    obs_a.dependencies = {f_bd};

    LambdaObservableConfig obs_b = LambdaObservableConfig::scalar(
        "STAT_DEV_OBS_B",
        [group, c_b, f_bs](LambdaDecay& ctx, ObservableId id) {
            (void)id;
            const double f = std::real(ctx.FLAVOR()(f_bs, DataType::VALUE));
            const double c = std::real(ctx.W().getFR(group, c_b, QCDOrder::LO, ContributionType::SM));
            return c - 0.5 * f;
        }
    );
    obs_b.aliases = {"stat-dev-obs-b"};
    obs_b.flha = LhaID(994001, 2);
    obs_b.dependencies = {f_bs};

    decay.observables = {obs_a, obs_b};
    oi->add_lambda_decay(decay, true);
    oi->enable_obs();

    ObservableId id_a = ObservableMapper::id_of("stat-dev-obs-a");
    ObservableId id_b = ObservableMapper::id_of("stat-dev-obs-b");

    print_param_set("ObservableInterface dependencies for STAT_DEV_OBS_A:", oi->get_all_ops_deps(id_a));
    print_param_set("ObservableInterface dependencies for STAT_DEV_OBS_B:", oi->get_all_ops_deps(id_b));

    // This is the exact adapter path used by StatisticManager: it sees an IModel,
    // not the concrete ObservableInterface. If this prints the dependencies, the
    // Statistic layer can discover them too.
    auto stat_param_optimizer = std::make_shared<StatParamOptimizerProxy>();
    ObservableInterfaceProxy model(oi, stat_param_optimizer);

    print_param_set("IModel/Statistic adapter dependencies for STAT_DEV_OBS_A:", model.get_obs_deps(id_a));
    print_param_set("IModel/Statistic adapter dependencies for STAT_DEV_OBS_B:", model.get_obs_deps(id_b));

    print_predictions("Predictions at the current point:", oi->compute_all());

    // Mutate the three declared inputs through the same optimized path used in fits/scans.
    // This proves that the custom observables respond to StatParamOptimizerProxy updates.
    std::map<ParamId, double> shifted_inputs = {
        {f_bd, 0.190},
        {f_bs, 0.235},
        {m_top, 173.2}
    };
    print_predictions("Predictions after a Statistic-style parameter update:",
                      model.predict_optimized({}, shifted_inputs));

    // Full StatisticInterface construction also works. This does not run a fit because the
    // toy custom observables have no experimental FOBS entries by default, but it exercises
    // the same dependency discovery used by StatisticManager::update_cache().
    StatisticConfig stat_config;
    stat_config.MC_draws = 20;
    StatisticInterface stat(stat_config, oi);

    print_param_map("StatisticManager active nuisance/input parameters:",
                    stat.get_active_observable_dependencies());

    return 0;
}
