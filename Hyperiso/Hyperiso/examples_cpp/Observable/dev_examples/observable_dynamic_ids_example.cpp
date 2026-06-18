#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ObservableInterface.h"
#include "mapper_hub.hpp"

namespace {

void print_values(const std::vector<ObservableValue>& values) {
    for (const auto& value : values) {
        std::cout << value.id.str() << " = " << value.value;
        if (value.bin.has_value()) {
            std::cout << " in [" << value.bin->first << ", " << value.bin->second << "]";
        }
        std::cout << "\n";
    }
}

// Helper dev: add builtin observables using runtime names instead of hard-coding Observables enums.
// For builtin observables, ObservableMapper::id_of(name) gives an ObservableId that can be
// passed directly to ObservableInterface::add_observable(ObservableId, ...).
void add_observables_by_name(ObservableInterface& interface,
                             const std::vector<std::pair<std::string, QCDOrder>>& observables,
                             bool add_dependencies = true)
{
    for (const auto& [name, order] : observables) {
        ObservableId id = ObservableMapper::id_of(name);
        std::cout << "Adding builtin observable from dynamic id: " << ObservableMapper::str(id) << "\n";
        interface.add_observable(id, order, add_dependencies);
    }
}

} // namespace

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    HyperisoConfig config;
    config.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    // Make sure the dynamic mapper layer knows all builtin names/aliases.
    init_all_builtins();

    ObservableInterface interface;

    // 1) Builtin observables dynamically resolved by name.
    add_observables_by_name(interface, {
        {"BR_Bs__mu_mu", QCDOrder::NNLO},
        {"BR_Bd__mu_mu", QCDOrder::NNLO},
        {"BR_B__Xs_gamma", QCDOrder::NNLO}
    });

    // Binned observable using ObservableId instead of the static enum.
    ObservableId fl_id = ObservableMapper::id_of("F_L_B0__K*0_mu_mu");
    BinnedObservableId fl_bin(fl_id, {1.1, 6.0});
    interface.add_observable(fl_bin, QCDOrder::NNLO, true);

    interface.enable_obs();

    std::cout << "\nBuiltin observables added through dynamic ids:\n";
    print_values(interface.compute_observable(ObservableMapper::id_of("BR_Bs__mu_mu")));
    print_values(interface.compute_observable(ObservableMapper::id_of("BR_Bd__mu_mu")));
    print_values({interface.compute_observable(fl_bin)});

    // 2) True custom observable + true custom decay.
    // A custom observable must always belong to a decay. The lambda-decay API registers
    // the DecayId and ObservableId and routes computation to the lambda below.
    LambdaDecayConfig custom_decay;
    custom_decay.canonical = "DEV_DYNAMIC_DECAY";
    custom_decay.aliases = {"dev-dynamic-decay"};
    custom_decay.matching_scale = 81.0;
    custom_decay.hadronic_scale = 4.8;
    custom_decay.order = QCDOrder::LO;

    LambdaObservableConfig custom_obs = LambdaObservableConfig::scalar(
        "DEV_DYNAMIC_OBS",
        [](LambdaDecay& ctx, ObservableId id) {
            (void)ctx;
            std::cout << "Computing custom dynamic observable: " << id.str() << "\n";
            return 42.0;
        }
    );
    custom_obs.aliases = {"dev-dynamic-obs"};
    custom_obs.flha = LhaID(990200, 1);
    custom_decay.observables.push_back(custom_obs);

    interface.add_lambda_decay(custom_decay, true);

    ObservableId custom_id = ObservableMapper::id_of("dev-dynamic-obs");
    auto parent = DecayMapper::get_decay_id(custom_id);
    if (parent) {
        std::cout << "\nCustom observable is linked to dynamic decay: " << parent->str() << "\n";
    }

    std::cout << "Custom observable value:\n";
    print_values(interface.compute_observable(custom_id));

    return 0;
}
