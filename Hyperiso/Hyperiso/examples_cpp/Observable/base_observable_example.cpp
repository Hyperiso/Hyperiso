#include <iostream>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ObservableInterface.h"
#include "BKsllDecay.h"
#include "BVFFCalculator.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    HyperisoConfig config;
    config.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    // Creation of the Observable interface. This interface needs to be created after the HyperisoMaster.
    // It will register all observables needed for the calculation.
    ObservableInterface interface;

    // For some decays, you can fine-tune the config if you want to use a specific set of form factors or change calculation parameters.
    // For heavy calculations (e.g. Kstarll) with integration, you can use kstar_conf.n_threads to run on multiple threads.
    // This is also available within the interface with interface.set_bkstarll_threads(nb_of_threads), same for bkll and bsphi.
    BKstarllConfig kstar_conf;

    // Use add_observable to add a new observable, with its order in QCD for the calculation.
    interface.add_observable(Observables::BR_B_XS_GAMMA, QCDOrder::NNLO);

    // With bin, use the BinnedObservableId class to add the bin. You can add multiple bins, which will be treated as different observables.
    interface.add_observable(BinnedObservableId(Observables::F_L_B0__KSTAR0_MU_MU, {1.1, 6.0}), QCDOrder::NNLO);

    // If you want to add all the observables of a decay, use the add_observables(Decays, QCDOrder) API.
    interface.add_observables(Decays::B__l_l, QCDOrder::NNLO);

    // Compute the observable. The API returns a vector of ObservableValue, because an observable can have several bins.
    for (const auto& value : interface.compute_observable(Observables::F_L_B0__KSTAR0_MU_MU)) {
        std::cout << value.id.str() << " = " << value.value;
        if (value.bin.has_value()) {
            std::cout << " in [" << value.bin->first << ", " << value.bin->second << "]";
        }
        std::cout << "\n";
    }

    // This allows to change the form factor choice of a particular decay.
    kstar_conf.ff_src = BV_FF_Src::GRvDV;

    // Set the new config inside the interface.
    interface.set_decay_config(Decays::B__Kstar_l_l, kstar_conf);

    std::cout << "\nAfter changing B->K* form-factor source:\n";
    for (const auto& value : interface.compute_observable(Observables::F_L_B0__KSTAR0_MU_MU)) {
        std::cout << value.id.str() << " = " << value.value << "\n";
    }

    for (const auto& value : interface.compute_observable(Observables::BR_BS_MUMU)) {
        std::cout << value.id.str() << " = " << value.value << "\n";
    }

    // This API allows to calculate all observables present in the interface.
    std::cout << "\nAll registered observables:\n";
    for (const auto& [obs, values] : interface.compute_all()) {
        std::cout << obs.str() << "\n";
        for (const auto& value : values) {
            std::cout << "  " << value.value << "\n";
        }
    }

    return 0;
}
