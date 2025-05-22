#include "HyperisoMaster.h"
#include "WilsonBuilder.h"
#include <iostream>
#include <cassert>

int main(){
    HyperisoMaster hyp = HyperisoMaster();
    Config config;
    config.model = Model::SUSY;

    hyp.init("lha/testInput.slha", config);

    WilsonBuilder builder;
    WilsonBuildConfig wilson_config;
    wilson_config.groups = {WGroup::B, WGroup::BPrime};
    wilson_config.matching_scale = 85.0;
    wilson_config.hadronic_scale = 4.5;
    wilson_config.order = QCDOrder::NLO;
    builder.build(wilson_config);

    // BlockProxy().log_all_blocks(ParameterType::WILSON);

    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD");
    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_TRADITIONAL");
    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_EW_SCALE");

    WilsonProvider wilson_provider = *builder.get_wilson_provider();

    std::shared_ptr<WilsonRequest> C7_request = std::make_shared<WilsonRequest>(
        WGroup::B,
        WCoef::C7,
        QCDOrder::NLO,
        ContributionType::SM,
        ScaleType::MATCHING,
        false
    );

    scalar_t C7w_SM = wilson_provider.get(C7_request);
    C7_request->contribution = ContributionType::BSM;
    scalar_t C7w_BSM = wilson_provider.get(C7_request);
    C7_request->contribution = ContributionType::TOTAL;
    scalar_t C7w_TOTAL = wilson_provider.get(C7_request);
    
    C7_request->scale_type = ScaleType::HADRONIC;
    C7_request->contribution = ContributionType::SM;
    scalar_t C7b_SM = wilson_provider.get(C7_request);
    C7_request->contribution = ContributionType::BSM;
    scalar_t C7b_BSM = wilson_provider.get(C7_request);
    C7_request->contribution = ContributionType::TOTAL;
    scalar_t C7b_TOTAL = wilson_provider.get(C7_request);

    C7_request->sum_qcd_orders = true;
    C7_request->contribution = ContributionType::SM;
    scalar_t C7b_SM_full = wilson_provider.get(C7_request);
    C7_request->contribution = ContributionType::BSM;
    scalar_t C7b_BSM_full = wilson_provider.get(C7_request);
    C7_request->contribution = ContributionType::TOTAL;
    scalar_t C7b_TOTAL_full = wilson_provider.get(C7_request);

    LOG_INFO("C7_MATCHING_SM =", C7w_SM);
    LOG_INFO("C7_MATCHING_BSM =", C7w_BSM);
    LOG_INFO("C7_MATCHING_TOTAL =", C7w_TOTAL);
    LOG_INFO("C7_HADRONIC_SM =", C7b_SM);
    LOG_INFO("C7_HADRONIC_BSM =", C7b_BSM);
    LOG_INFO("C7_HADRONIC_TOTAL =", C7b_TOTAL);
    LOG_INFO("C7_HADRONIC_FULL_SM =", C7b_SM_full);
    LOG_INFO("C7_HADRONIC_FULL_BSM =", C7b_BSM_full);
    LOG_INFO("C7_HADRONIC_FULL_TOTAL =", C7b_TOTAL_full);

    assert(std::abs(C7w_SM + C7w_BSM - C7w_TOTAL) < 1e-15);
    assert(std::abs(C7b_SM + C7b_BSM - C7b_TOTAL) < 1e-15);

    return 0;
}
