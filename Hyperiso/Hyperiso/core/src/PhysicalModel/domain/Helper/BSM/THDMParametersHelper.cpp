#include "THDMParametersHelper.h"


void thdm_parameters::init(int gen, WGroupId grp) {
    if (initialized) {
		return;
	}

    init_scale_independent_block(gen);
    init_matching_block();

    initialized = true;
}

void thdm_parameters::init_scale_independent_block(int gen) {
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
        int id {1};
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, xh, 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, xH, 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, xH0, 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, xA, 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, m_H, 0., 0.));
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, beta, 0., 0.)); 
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, lu, 0., 0.));
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, ld, 0., 0.));
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, alpha, 0., 0.));
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, le, 0., 0.)); 
    };

    iblock_c->compose_block("WPARAM_SI_BSM", src, func);
 }

void thdm_parameters::init_matching_block() {

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