#include "Parameters.h"
#include "HyperisoMaster.h"
#include "BlockProvider.h"
#include "CorrelationProvider.h"
#include <iostream>
#include <cmath>

int main() {
    std::cout << "here ! " << std::endl;
    auto hyp = HyperisoMaster();  // Initialize program manager with LHA file containing SMINPUTS block
    HyperisoConfig config;
    config.model = Model::SM;
    std::cout << "here ! " << std::endl;
    hyp.init("Test/InputFiles/testInput.slha", config);
    auto sm_params = Parameters::GetInstance(ParameterType::SM); // SM Model
    // auto susy_params = Parameters::GetInstance(ParameterType::BSM); // SUSY Model
    
    std::shared_ptr<BlockAccessor> truc = sm_params->get_block_accessor();

    auto mmh =  truc->get_block_sources("VCKM");

    for (auto mm : mmh) {
        std::cout << mm.first << std::endl;
    }
    auto ooooh = truc->get_all_source_parameters({ParamId(ParameterType::SM, "VCKM", {2,2})});

    for (auto oo : ooooh) {
        std::cout << oo << std::endl;
    }
    try {
        
        double alpha_s_MZ = pow((*sm_params)("GAUGE", 1), 2) / (4 * M_PI);
        std::cout << "Alpha_s(MZ): " << alpha_s_MZ << std::endl;

        std::cout << "mass up quark : " << (*sm_params)("MASS", 2) << std::endl;
        std::cout << "first element of ckm (real part) : " << (*sm_params)("RECKM", 0) << std::endl;
        //  double susy_mass = (*susy_params)("MASS", 1000021); // Example PDG code for a SUSY particle
        // std::cout << "SUSY Particle Mass: " << susy_mass << std::endl;


    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    // try {
    //     std::cout << (*susy_params)(std::string("GAUGE"), 1) << std::endl;
    // } catch(const std::exception& e) {
    //     std::cerr << "Error: " << e.what() << std::endl;
    // }
 
    // auto mm2 = MemoryManager::GetInstance("Test/InputFiles/testinput_thdm.lha", {0,2});
    // mm2->init();
    // auto thdm_params = Parameters::GetInstance(2); // THDM Model
    // std::cout << "THDM matrice yu : " << (*thdm_params)("YU", 22) << std::endl;
    BlockProvider().log_all_blocks(ParameterType::OBSERVABLE);
    LOG_INFO(CorrelationProvider()("DEFAULT", Observables::BR_BS_MUMU_UNTAG, Observables::BR_BD_MUMU, CorrelationProvider::CorrelationType::COMBINED));
    LOG_INFO(CorrelationProvider()("DEFAULT2", Observables::BR_BS_MUMU_UNTAG, Observables::BR_BD_MUMU, CorrelationProvider::CorrelationType::COMBINED));
    return 0;
}