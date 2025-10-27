#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"
#include "BlockProxy.h"

int main() {
    auto hyp = HyperisoMaster();
    Config config_hyp;
    hyp.init("Test/InputFiles/testInput.flha", config_hyp); // Initialize program manager with LHA file

    auto wi = WilsonInterface(); // Initialize interface and build the required groups

    WilsonBuildConfig config({WGroup::B}, 81, 4.18, QCDOrder::NNLO);

    wi.build(config);

    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD");
    config.groups = {GroupMapper::to_id(WGroup::BPrime)};
    wi.addWilsonGroup(config);

    BlockProxy().log_block(ParameterType::WILSON, "BPrimeCoefficients_B_SCALE_STANDARD");


    return 0;
}