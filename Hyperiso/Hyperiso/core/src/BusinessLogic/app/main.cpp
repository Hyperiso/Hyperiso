#include "Logger.h"
#include <iostream>
#include <cassert>
#include "ObservableInterface.h"
#include "HyperisoMaster.h"

// Utility function for displaying results
template<typename T, typename U>
void assert_equal(const T& result, const U& expected, const std::string& test_name) {
    if (result == expected) {
        std::cout << "[PASS] " << test_name << "\n";
    } else {
        std::cerr << "[FAIL] " << test_name << ": got " << result << ", expected " << expected << "\n";
    }
}

void test_interface(ObservableInterface& interface) {
    using enum Observables;
    using enum QCDOrder;
    using enum UncertaintyType;

    // Add single observable
    interface.add_observable(BR_BS_MUMU, NLO, true);
    assert_equal(interface.get_current_observables().count(ObservableMapper::to_id(BR_BS_MUMU)), 1u, "add_observable");

    // Add multiple observables (map version)
    interface.add_observables({
        {R_D, LO},
        {R_DSTAR, LO}
    }, true);
    assert_equal(interface.get_current_observables().count(ObservableMapper::to_id(R_D)), 1u, "add_observables (map) - R_D");
    assert_equal(interface.get_current_observables().count(ObservableMapper::to_id(R_DSTAR)), 1u, "add_observables (map) - R_DSTAR");

    // Add decay observables
    interface.add_observables(Decays::B__l_l, LO);
    for (auto obs : DecayMapper::get_observables(Decays::B__l_l)) {
        assert_equal(interface.get_current_observables().count(ObservableMapper::to_id(obs)), 1u, "add_observables (decay)");
    }

    // Add observable parameters
    // interface.add_observable_parameter(R_DSTAR, ParamId(ParameterType::FLAVOR, "FMASS", 423));
    // interface.add_observable_parameters(R_D, {ParamId(ParameterType::FLAVOR, "FMASS", 421), ParamId(ParameterType::FLAVOR, "FMASS", 521)});

    // Compute observable
    auto val = interface.compute_observable(R_DSTAR);
    for (auto elem : val) {
        std::cout << "compute_observable: " << elem.value << "\n";
    }

    // Compute uncertainty
    auto unc = interface.compute_uncertainty(R_DSTAR);
    std::cout << "compute_uncertainty: " << unc << "\n";

    // Leading uncertainties
    auto leading_uncs = interface.compute_leading_uncertainties(BR_BS_MUMU, 2);
    std::cout << "compute_leading_uncertainties size: " << leading_uncs.size() << "\n";

    // All uncertainties
    auto all_uncs = interface.compute_all_uncertainties();
    std::cout << "compute_all_uncertainties size: " << all_uncs.size() << "\n";

    // Chi2
    double chi2 = interface.compute_chi2();
    std::cout << "compute_chi2: " << chi2 << "\n";

    // Remove single observable
    interface.remove_observable(BR_BS_MUMU);
    assert_equal(interface.get_current_observables().count(ObservableMapper::to_id(BR_BS_MUMU)), 0u, "remove_observable");

    // Remove multiple observables (set)
    interface.remove_observables({R_D, R_DSTAR});
    assert_equal(interface.get_current_observables().count(ObservableMapper::to_id(R_D)), 0u, "remove_observables (set) - R_D");

    // Remove decay observables
    interface.remove_observables(Decays::B__l_l);
    for (auto obs : DecayMapper::get_observables(Decays::B__l_l)) {
        assert_equal(interface.get_current_observables().count(ObservableMapper::to_id(obs)), 0u, "remove_observables (decay)");
    }

    // Experimental values and uncertainties
    interface.add_observable(R_DSTAR, LO); // Add back for exp tests
    scalar_t exp_val = interface.get_exp_value(R_DSTAR);
    scalar_t stat = interface.get_exp_uncertainty(R_DSTAR, STAT);
    scalar_t syst = interface.get_exp_uncertainty(R_DSTAR, SYST);
    std::cout << "get_exp_value: " << exp_val << ", stat: " << stat << ", syst: " << syst << "\n";

    // compute_all
    auto estimates = interface.compute_all();
    std::cout << "compute_all size: " << estimates.size() << "\n";

    // get_all_exp
    auto all_exp = interface.get_all_exp();
    std::cout << "get_all_exp size: " << all_exp.size() << "\n";
}


int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    HyperisoMaster hyperiso;
    Config config;
    config.model = Model::SM;
    hyperiso.init("default/lha/testInput.flha", config);
    LOG_INFO("HyperisoMaster initialized");

    auto interface = ObservableInterface(); 

    test_interface(interface);
}