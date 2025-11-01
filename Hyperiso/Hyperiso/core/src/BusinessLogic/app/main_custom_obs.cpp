#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"
#include "Wilson/WilsonObs.h"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyp;
    Config config;
    config.model = Model::SM;
    hyp.init("lha/si_input.flha", config);
    LOG_INFO("HyperisoMaster initialized");

    QCDOrder order = QCDOrder::NNLO;
    ObservableInterface oi;

    DecayMapper::register_custom("Wilson");
    auto wilson_decay = std::make_shared<WilsonDecay>(DecayMapper::id_of("Wilson"), order, 81, 5);
    oi.add_custom_decay(DecayMapper::id_of("Wilson"), wilson_decay);
    ObservableId o = ObservableMapper::id_of("C7");
    oi.add_observable(o, order, false);
    oi.add_observable_parameters(o, {
        ParamId {ParameterType::SM, "SMINPUTS", 3},
        ParamId {ParameterType::SM, "SMINPUTS", 4},
        ParamId {ParameterType::SM, "SMINPUTS", 6},
        ParamId {ParameterType::SM, "MASS", 24}
    });

    LOG_INFO("C7 =", oi.compute_observable(o)[0].value, "+-", oi.compute_uncertainty(o));

    auto print_leading = [&oi] (ObservableId o, size_t n) {
        LOG_INFO("---------- Leading uncertainties for", ObservableMapper::str(o));
        for (const auto& [pid, u] : oi.compute_leading_uncertainties(o, n)) {
            LOG_INFO("\t-", pid, ":", u.real());
        }
    };

    print_leading(o, 4);

    return 0;
}