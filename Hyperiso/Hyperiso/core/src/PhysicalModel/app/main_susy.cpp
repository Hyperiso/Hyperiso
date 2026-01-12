#include <iostream>

#include "MemoryManager.h"
#include "Parameters.h"
#include "WilsonInterface.h"

int main() {
    auto hyp = HyperisoMaster();
    HyperisoConfig config_hyp;
    // config_hyp.flags[ExternalFlag::USE_MARTY] = false;
    config_hyp.model = Model::SUSY;
    config_hyp.mty_model_name = "ZPrime";
    config_hyp.mty_model_path = "/home/theo/hyperiso/Assets/input_files/marty_model/ZPrime.h";

    hyp.init("/home/theo/hyperiso/Hyperiso/Hyperiso/core/Test/InputFiles/testInput.slha", config_hyp); // Initialize program manager with LHA file

    auto wi = WilsonInterface(); // Initialize interface and build the required groups

    WilsonBuildConfig config({WGroup::BScalar, WGroup::B, WGroup::BPrime}, 81, 42, QCDOrder::NNLO);

    wi.build(config);

    // BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD");

    for (auto elem : {WCoef::CQ1, WCoef::CQ2}) {
        LOG_INFO(WCoefMapper::str(elem), "(MW) at LO =", wi.getM(WGroup::BScalar, elem, QCDOrder::LO, ContributionType::BSM));
        LOG_INFO(WCoefMapper::str(elem), "(MW) at NLO =", wi.getM(WGroup::BScalar, elem, QCDOrder::NLO, ContributionType::BSM));
        LOG_INFO("\n");
    }

    for (auto elem : {WCoef::C1, WCoef::C2, WCoef::C3, WCoef::C4, WCoef::C5, WCoef::C6, WCoef::C7, WCoef::C8, WCoef::C9, WCoef::C10}) {
        LOG_INFO(WCoefMapper::str(elem), "(MW) at LO =", wi.getM(WGroup::B, elem, QCDOrder::LO, ContributionType::BSM));
        LOG_INFO(WCoefMapper::str(elem), "(MW) at NLO =", wi.getM(WGroup::B, elem, QCDOrder::NLO, ContributionType::BSM));
        LOG_INFO(WCoefMapper::str(elem), "(MW) at NNLO =", wi.getM(WGroup::B, elem, QCDOrder::NNLO, ContributionType::BSM));
        LOG_INFO("\n");
    }

    for (auto elem : {WCoef::CP1, WCoef::CP2, WCoef::CP3, WCoef::CP4, WCoef::CP5, WCoef::CP6, WCoef::CP7, WCoef::CP8, WCoef::CP9, WCoef::CP10, WCoef::CPQ1, WCoef::CPQ2}) {
        LOG_INFO(WCoefMapper::str(elem), "(MW) at LO =", wi.getM(WGroup::BPrime, elem, QCDOrder::LO, ContributionType::BSM));
        LOG_INFO(WCoefMapper::str(elem), "(MW) at NLO =", wi.getM(WGroup::BPrime, elem, QCDOrder::NLO, ContributionType::BSM));
        LOG_INFO("\n");
    }

    LOG_INFO("Parameters created");
    // LOG_INFO("C9(mu_h) at LO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::LO, ContributionType::TOTAL));
    // LOG_INFO("AGAIN");
    // LOG_INFO("C9(mu_h) at NLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NLO, ContributionType::TOTAL));
    // LOG_INFO("C9(mu_h) at NNLO =", wi.getR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    
    // LOG_INFO("C9(mu_h) full =", wi.getFR(WGroup::B, WCoef::C9, QCDOrder::NNLO, ContributionType::TOTAL));
    return 0;
}