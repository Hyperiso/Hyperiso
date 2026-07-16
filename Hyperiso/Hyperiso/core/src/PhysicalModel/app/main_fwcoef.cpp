#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"
#include "BlockProxy.h"

int main() {
    auto hyp = HyperisoMaster();
    HyperisoConfig config_hyp;
    config_hyp.flags[ExternalFlag::HAS_WILSON_INPUT] = true;
    hyp.init("Test/InputFiles/testInput.flha", config_hyp); // Initialize program manager with LHA file

    auto wi = WilsonInterface(); // Initialize interface and build the required groups

    WilsonBuildConfig config({WGroup::B}, 81, 4.18, QCDOrder::NLO);

    wi.build(config);
    BlockProxy().log_block(ParameterType::WILSON, "FWCOEF");
    auto list = BlockProxy().get_block_list(ParameterType::WILSON);
    for (auto elem : list) {
        std::cout << elem << std::endl;
    }
    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_EW_SCALE");
    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_EW_SCALE_SM");
    BlockProxy().log_all_blocks(ParameterType::WILSON);
    std::cout << "here !" << std::endl;

    std::vector<std::string> name {"C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10"};
    for (auto& coeff : name) {
        complex_t C = {0,0};
            C = wi.getMatchingCoefficient(WGroup::B,WCoefMapper::enum_of(WCoefMapper::enum_elt(coeff)).value(), QCDOrder::LO, ContributionType::SM );

        std::cout << "," << C.real() << "," << C.imag() << std::endl;
    }
    BlockProxy().log_all_blocks(ParameterType::WILSON);
    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD");



    return 0;
}