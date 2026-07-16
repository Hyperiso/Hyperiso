#include <iostream>
#include <vector>

#include "HyperisoMaster.h"
#include "Include.h"
#include "Logger.h"
#include "ObservableInterface.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);

    // Dev example: manually scan one input parameter and recompute an observable.
    HyperisoConfig config;
    config.model = Model::SM;

    HyperisoMaster hyp;
    hyp.init("lha/si_input.flha", config);

    ObservableInterface interface;
    interface.add_observable(Observables::BR_BS_MUMU, QCDOrder::NNLO, true);

    const std::vector<double> fbs_values = {0.200, 0.225, 0.250};

    for (double fbs : fbs_values) {
        // set_param is a convenience wrapper around the global Parameters storage.
        // Here FCONST[531,1] is used as an example fit/scan parameter.
        interface.set_param("FCONST", LhaID(531, 1), fbs, ParameterType::FLAVOR);

        // After changing parameters, reload the cached decay parameters before computing again.
        interface.reload_params();
        interface.enable_obs();

        const auto values = interface.compute_observable(Observables::BR_BS_MUMU);
        std::cout << "FCONST[531,1] = " << fbs << " -> BR_BS_MUMU = " << values.front().value << "\n";
    }

    return 0;
}
