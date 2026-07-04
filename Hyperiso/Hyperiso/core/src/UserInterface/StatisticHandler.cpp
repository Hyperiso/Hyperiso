#include "StatisticHandler.h"

#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "CliUtils.h"
#include "ObservableInterface.h"
#include "StatisticInterface.h"
#include "mapper_hub.hpp"

namespace {

void print_statistic_usage() {
    std::cout
        << "Usage:\n"
        << "  hyperiso-ui statistic summary [options]\n\n"
        << "Options:\n"
        << "  --observables <csv>     Observable names, default BR_Bs__mu_mu,BR_B__Xs_gamma\n"
        << "  --draws <n>             MC draws for uncertainty mode, default 200\n"
        << "  --seed <n>              RNG seed for reproducible MC runs, default 123456\n"
        << "  --uncertainties         Also compute MC uncertainty summaries\n"
        << "  --chi2                  Use CHI2_MC_COVARIANCE advanced likelihood mode\n"
        << "  --progress              Show MC progress bar when MC is run\n"
        << "  --samples-csv <path>    Write accepted MC observable samples to CSV\n"
        << "  --order <order>         LO, NLO or NNLO, default NNLO\n"
        << "  --model <model>         SM, THDM, MSSM or MARTY, default SM\n"
        << "  --lha <path>            Input LHA/FLHA file\n";
}

void print_predictions(const std::map<ObservableId, std::vector<ObservableValue>>& predictions) {
    for (const auto& [obs, values] : predictions) {
        std::cout << "\n" << ObservableMapper::str(obs) << "\n";
        for (const auto& value : values) {
            std::cout << "  prediction=" << value.value;
            if (value.bin.has_value()) {
                std::cout << " in [" << value.bin->first << ", " << value.bin->second << "]";
            }
            std::cout << "\n";
        }
    }
}

} // namespace

int handleStatisticOptions(int argc, char* argv[]) {
    CliOptions opts = CliOptions::parse(argc, argv, 1);
    const std::string command = opts.positionals.empty() ? "summary" : opts.positionals[0];

    if (opts.flag("help", false) || command == "help") {
        print_statistic_usage();
        return 0;
    }
    if (command != "summary") {
        throw std::invalid_argument("Unknown statistic command: " + command);
    }

    auto hyp = init_hyperiso_from_cli(opts);
    init_all_builtins();

    const QCDOrder order = parse_qcd_order(opts.get("order", "NNLO"));
    const auto obs_names = opts.list("observables", {"BR_Bs__mu_mu", "BR_B__Xs_gamma"});

    auto oi = std::make_shared<ObservableInterface>();
    for (const auto& name : obs_names) {
        oi->add_observable(ObservableMapper::id_of(name), order, true);
    }
    oi->enable_obs();

    StatisticConfig cfg;
    cfg.MC_draws = static_cast<std::size_t>(opts.get_int("draws", 200));
    cfg.MC_seed = static_cast<unsigned int>(opts.get_int("seed", 123456));
    cfg.print_mc_progress = opts.flag("progress", false);
    if (opts.has("samples-csv")) {
        cfg.write_mc_samples_csv = true;
        cfg.mc_samples_csv_path = opts.get("samples-csv");
    }
    cfg.print_fit_summary = opts.flag("verbose", false);
    cfg.print_scan_summary = opts.flag("verbose", false);
    cfg.advanced.likelihood_mode = opts.flag("chi2", false)
        ? StatisticLikelihoodMode::CHI2_MC_COVARIANCE
        : StatisticLikelihoodMode::PROFILED_NUISANCE;

    StatisticInterface stat(cfg, oi);

    print_section("Statistic summary");
    std::cout << "observables=" << obs_names.size()
              << ", draws=" << cfg.MC_draws
              << ", seed=" << cfg.MC_seed
              << ", mode=" << (opts.flag("chi2", false) ? "CHI2_MC_COVARIANCE" : "PROFILED_NUISANCE")
              << "\n";

    print_section("Predictions");
    print_predictions(oi->compute_all());

    print_section("Active observable dependencies");
    const auto deps = stat.get_active_observable_dependencies();
    if (deps.empty()) {
        std::cout << "  <none>\n";
    } else {
        for (const auto& [pid, value] : deps) {
            std::cout << "  " << pid << " = " << value << "\n";
        }
    }

    if (opts.flag("uncertainties", false)) {
        print_section("MC uncertainty summaries");
        const auto summaries = stat.compute_uncertainties();
        for (const auto& [obs, summary] : summaries) {
            std::cout << "  " << obs.str() << " -> " << summary << "\n";
        }
    }

    return 0;
}
