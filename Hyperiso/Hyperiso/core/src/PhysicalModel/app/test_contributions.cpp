#include "HyperisoMaster.h"
#include "WilsonBuilder.h"
#include <iostream>
#include <cassert>

int main() {
    HyperisoMaster hyp = HyperisoMaster();
    Config config;
    config.model = Model::CUSTOM;
    config.flags[ExternalFlag::USE_MARTY] = true; // TODO : Théo not happy
    config.mty_model_name = "THDM_Model";
    config.mty_model_path = std::string(project_tp_root.data()) + "MARTY/src/MARTY/src/marty/models/thdm.h";

    hyp.init("lha/testinput_thdm.lha", config);

    WilsonBuildConfig wilson_config;
    // wilson_config.groups = {WGroup::B, WGroup::BPrime};
    wilson_config.groups = {WGroup::B};
    wilson_config.matching_scale = 85.0;
    wilson_config.hadronic_scale = 4.5;
    wilson_config.order = QCDOrder::NNLO;
    WilsonBuilder builder {wilson_config};
    // builder.build(wilson_config);

    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_STANDARD");
    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_B_SCALE_TRADITIONAL");
    BlockProxy().log_block(ParameterType::WILSON, "BCoefficients_EW_SCALE");

    WilsonProvider wilson_provider = *builder.get_wilson_provider();

    std::shared_ptr<WilsonRequest> C7_request = std::make_shared<WilsonRequest>(
        WGroup::B,
        WCoef::C5,
        QCDOrder::NNLO,
        ContributionType::SM,
        ScaleType::MATCHING,
        false
    );
    C7_request->basis = WilsonBasis::B_STANDARD;

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

    assert(std::abs(C7w_SM + C7w_BSM - C7w_TOTAL) < 1e-13);
    assert(std::abs(C7b_SM + C7b_BSM - C7b_TOTAL) < 1e-13);

    return 0;
}
