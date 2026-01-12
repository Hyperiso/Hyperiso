#include "HyperisoMaster.h"
#include "WilsonBuilder.h"
#include <iostream>
#include <cassert>

int main() {
    HyperisoMaster hyp = HyperisoMaster();
    HyperisoConfig config;
    config.model = Model::THDM;
    // config.flags[ExternalFlag::USE_MARTY] = false; // TODO : Théo not happy
    config.mty_model_name = "THDM";
    config.mty_model_path = std::string(project_tp_root.data()) + "MARTY/src/MARTY/src/marty/models/thdm.h";

    hyp.init("lha/testinput_thdm.lha", config);

    WilsonBuildConfig wilson_config;
    wilson_config.groups = {GroupMapper::to_id(WGroup::B), GroupMapper::to_id(WGroup::BPrime)};
    wilson_config.matching_scale = 80;
    wilson_config.hadronic_scale = ParameterProvider()({ParameterType::SM, "QCD", {5, 3}}) / 2;
    wilson_config.order = QCDOrder::LO;
    WilsonBuilder builder {wilson_config};

    

    return 0;
}
