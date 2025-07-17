#include "Logger.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <filesystem>
#include <cassert>

#include "ObservableInterface.h"
#include "HyperisoMaster.h"
#include "config.hpp"

int main() {
    Logger::getInstance()->setLevel(Logger::LogLevel::INFO);
    namespace fs = std::filesystem;

    std::string lha_dir = project_assets_root.data() + std::string("scan");
    std::ofstream ofs("B_s_gamma_SI_scan_LO.csv");
    ofs << "lha_file,M_Hp,BR,u(BR)\n";

    // Prépare HyperIso
    Config config;
    config.model = Model::THDM;
    config.flags[ExternalFlag::USE_MARTY] = false;
    config.mty_model_name = "THDM";
    config.mty_model_path = project_assets_root.data() + std::string("input_files/marty_model/thdm.h");
    HyperisoMaster hyp;
    hyp.init("lha/testinput_thdm.lha", config);
    
    // Parcours tous les fichiers LHA
    for (const auto& entry : fs::directory_iterator(lha_dir)) {
        if (entry.path().extension() != ".lha") continue;

        std::string lha_path = entry.path().string();
        std::string lha_name = entry.path().filename().string();

        LOG_INFO("Processing " + lha_name);

        
        hyp.init(lha_path, config);
        MemoryManager::GetInstance()->switch_lha(lha_path, config);
        auto sm = SMModelStrategy();
        std::cout << "aah" << std::endl;
        sm.add_absent_block(SMModelStrategy().initializeParameters(*Parameters::GetInstance(ParameterType::SM)));
        sm.postInitialization(*Parameters::GetInstance(ParameterType::SM));
        BSMModelStrategy().initializeParameters(*Parameters::GetInstance(ParameterType::BSM));
        // Parameters::GetInstance(ParameterType::SM)->init_blocks(ParameterType::SM);
        // SMModelStrategy().postInitialization(*Parameters::GetInstance(ParameterType::SM));
        // Parameters::GetInstance(ParameterType::BSM)->init_blocks(ParameterType::BSM);
        std::cout << "zh" << std::endl;
        ParameterProvider pp;
        // ParameterSetter ps;
        std::cout << "ah" << std::endl;
        // WilsonInterface wi;
        std::cout << "ch" << std::endl;
        WilsonBuildConfig wilson_config;
        std::cout << "dh" << std::endl;
        // wilson_config.groups = {WGroup::B, WGroup::BPrime};
        // wilson_config.matching_scale = 81;
        // wilson_config.hadronic_scale = pp({ParameterType::SM, "QCD", {5, 3}}) / 2;
        wilson_config.order = QCDOrder::LO;
        std::cout << "eh" << std::endl;
        // wi.build(wilson_config);
        std::cout << "hh" << std::endl;
        ObsWilsonHelper(true);
        ObservableInterface oi;
        std::cout << "bh" << std::endl;
        // ObservableInterface oi;
        oi.add_observable(Observables::BR_B_XS_GAMMA, wilson_config.order, true);

        // Extraire M_H+ depuis la LHA (paramètre avec PDG 37)
        double mhp = pp({ParameterType::BSM, "MASS", 37});
        // ps.mutate({ParameterType::BSM, "MASS", 37}, mhp);  // Assure la mise à jour

        double br = oi.compute_observable(Observables::BR_B_XS_GAMMA).real();
        double ubr = oi.compute_uncertainty(Observables::BR_B_XS_GAMMA).real();

        ofs << lha_name << "," << std::setprecision(8) << mhp << "," << br << "," << ubr << "\n";
        
    }

    return 0;
}