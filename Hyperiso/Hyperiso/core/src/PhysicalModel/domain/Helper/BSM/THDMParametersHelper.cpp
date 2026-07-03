#include "THDMParametersHelper.h"
#include "wcoef_ids.hpp"


void THDMParameterHelper::init(int gen, WGroupId) {
    if (initialized) {
		return;
	}

    init_scale_independent_block(gen);
    init_matching_block();

    initialized = true;
}

void THDMParameterHelper::init_scale_independent_block(int) {
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS"}}, {ParameterType::BSM, {"MASS", "ALPHA", "MINPAR", "YU", "YD", "YE"}},
                                                                        {ParameterType::WILSON, {"WPARAM_SI_SM"}}};

    auto func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
        double xh = pow(src.get_val("MASS", 25) / src.get_val("MASS", 24), 2);
        double alpha = src.get_val("ALPHA", LhaID(""));
        double m_H = src.get_val("MASS" ,37);
        double beta = atan(src.get_val("MINPAR" ,3));
        double lu = src.get_val("YU" ,LhaID(2,2));
        double ld = src.get_val("YD" ,LhaID(2,2));
        double gen = src.get_val("WPARAM_SI_SM" , 2);
        double le = src.get_val("YE" ,LhaID(gen-1, gen-1));
        double mW = src.get_val("MASS", 24);
        double xH=pow(m_H/mW,2.);
        double xH0=pow(src.get_val("MASS" ,35) / mW, 2.);
        double xA=pow(src.get_val("MASS" ,36) / mW, 2.);

        dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 1}, xh, 0., 0.));
		dep_block->store_or_assign(2, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 2}, xH, 0., 0.));
		dep_block->store_or_assign(3, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 3}, xH0, 0., 0.));
		dep_block->store_or_assign(4, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 4}, xA, 0., 0.));
		dep_block->store_or_assign(5, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 5}, m_H, 0., 0.));
        dep_block->store_or_assign(6, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 6}, beta, 0., 0.)); 
        dep_block->store_or_assign(7, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 7}, lu, 0., 0.));
        dep_block->store_or_assign(8, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 8}, ld, 0., 0.));
        dep_block->store_or_assign(9, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 9}, alpha, 0., 0.));
        dep_block->store_or_assign(10, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", 10}, le, 0., 0.));

        // Keep the historical slot 10 for the configured lepton generation, and
        // additionally expose all three lepton Yukawas for CQ/CPQ_E/MU/TA.
        for (int i = 0; i < 3; ++i) {
            const int slot = WCoefMapper::thdm_lepton_yukawa_slot_from_index(i);
            const double le_i = src.get_val("YE", LhaID(i, i));
            dep_block->store_or_assign(slot, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", slot}, le_i, 0., 0.));
        }
    };

    iblock_c->compose_block("WPARAM_SI_BSM", src, func);
 }

void THDMParameterHelper::init_matching_block() {

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::SM, {"MASS" /*, "QCD"*/}}, 
        {ParameterType::WILSON, {"WPARAM_MATCH_SM"}}
    };

    auto func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
        double yt = pow(src.get_val("WPARAM_MATCH_SM" ,6)/src.get_val("MASS" ,37),2.); // param->mass_H (25)
		dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}, yt, 0., 0.));
    };

    iblock_c->compose_block("WPARAM_MATCH_BSM", src, func);

}