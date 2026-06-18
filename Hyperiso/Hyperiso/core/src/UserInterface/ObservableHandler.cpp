#include "ObservableHandler.h"

#include <iostream>
#include <sstream>
#include <vector>

#include "CliUtils.h"
#include "ObservableInterface.h"
#include "mapper_hub.hpp"

namespace {

void print_observable_usage() {
    std::cout
        << "Usage:\n"
        << "  hyperiso-ui observable summary [options]\n\n"
        << "Options:\n"
        << "  --observables <csv>       Observable names, default BR_Bs__mu_mu,BR_B__Xs_gamma\n"
        << "  --bins <specs>            Binned specs OBS:min:max, comma-separated\n"
        << "  --order <order>           LO, NLO or NNLO, default NNLO\n"
        << "  --model <model>           SM, THDM, MSSM or MARTY, default SM\n"
        << "  --lha <path>              Input LHA/FLHA file\n";
}

struct BinSpec {
    std::string observable;
    double low = 0.0;
    double high = 0.0;
};

BinSpec parse_bin_spec(const std::string& raw) {
    std::stringstream ss(raw);
    std::string obs, lo, hi;
    if (!std::getline(ss, obs, ':') || !std::getline(ss, lo, ':') || !std::getline(ss, hi, ':')) {
        throw std::invalid_argument("Invalid bin spec '" + raw + "'. Expected OBS:min:max");
    }
    return {obs, std::stod(lo), std::stod(hi)};
}

void print_values(const std::vector<ObservableValue>& values) {
    for (const auto& value : values) {
        std::cout << "  " << value.id.str() << " = " << value.value;
        if (value.bin.has_value()) {
            std::cout << " in [" << value.bin->first << ", " << value.bin->second << "]";
        }
        std::cout << "\n";
    }
}

} // namespace

int handleObservableOptions(int argc, char* argv[]) {
    CliOptions opts = CliOptions::parse(argc, argv, 1);
    const std::string command = opts.positionals.empty() ? "summary" : opts.positionals[0];

    if (opts.flag("help", false) || command == "help") {
        print_observable_usage();
        return 0;
    }
    if (command != "summary") {
        throw std::invalid_argument("Unknown observable command: " + command);
    }

    auto hyp = init_hyperiso_from_cli(opts);
    init_all_builtins();

    const QCDOrder order = parse_qcd_order(opts.get("order", "NNLO"));
    const auto obs_names = opts.list("observables", {"BR_Bs__mu_mu", "BR_B__Xs_gamma"});
    const auto bin_specs = opts.list("bins", {});

    ObservableInterface interface;
    for (const auto& name : obs_names) {
        interface.add_observable(ObservableMapper::id_of(name), order, true);
    }
    for (const auto& raw : bin_specs) {
        const BinSpec spec = parse_bin_spec(raw);
        interface.add_observable(
            BinnedObservableId(ObservableMapper::id_of(spec.observable), {spec.low, spec.high}),
            order,
            true
        );
    }

    interface.enable_obs();

    print_section("Observable summary");
    std::cout << "order=" << OrderMapper::str(order) << "\n";
    for (const auto& [obs, values] : interface.compute_all()) {
        std::cout << "\n" << ObservableMapper::str(obs) << "\n";
        print_values(values);
    }

    return 0;
}
