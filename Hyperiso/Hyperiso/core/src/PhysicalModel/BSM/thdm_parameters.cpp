#include "thdm_parameters.h"


void thdm_parameters::init() {
    if (thdm_parameters::initialized) {
		return;
	}

	std::cout << "Initializing thdm_parameters" << std::endl;

    init_scale_independant_block();
    init_matching_block();
}

void thdm_parameters::init_scale_independant_block() {
    

    // ParameterProxy(ParameterType::SM);
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS"}}, {ParameterType::SM, {"MASS"}}, 
                                                                       {ParameterType::SM, {"ALPHA"}}, {ParameterType::SM, {"HMIX"}},
                                                                       {ParameterType::SM, {"YU"}}, {ParameterType::SM, {"YD"}}, {ParameterType::SM,  {"YL"}},
                                                                        {ParameterType::WILSON, {"WPARAM_SI_SM"}}};

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double xh = pow(src.at("MASS")->retrieve(25)->get_val() / src.at("MASS")->retrieve(24)->get_val(), 2);

        double alpha = src.at("ALPHA")->retrieve(0)->get_val();
        
        double m_H = src.at("MASS")->retrieve(37)->get_val();
        double beta = atan(src.at("HMIX")->retrieve(2)->get_val());

        double lu = src.at("YU")->retrieve(22)->get_val();
        double ld = src.at("YD")->retrieve(22)->get_val();
        double gen = src.at("WPARAM_SI_SM")->retrieve(2)->get_val();
        double le = src.at("YL")->retrieve(10*(gen-1)+gen-1)->get_val();
        double mW = src.at("MASS")->retrieve(24)->get_val();
        double xH=pow(m_H/mW,2.);
        double xH0=pow(src.at("MASS")->retrieve(35)->get_val() / mW, 2.);
        double xA=pow(src.at("MASS")->retrieve(36)->get_val() / mW, 2.);

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

    thdm_parameters::composer.compose_block("WPARAM_SI_BSM", src, func);

}

void thdm_parameters::init_matching_block() {

    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::SM, {"MASS" /*, "QCD"*/}}, 
        {ParameterType::WILSON, {"WPARAM_MATCH_SM"}}
    };

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double yt = pow(src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val()/src.at("MASS")->retrieve(37)->get_val(),2.); // param->mass_H (25)
		dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}, yt, 0., 0.));

		LOG_INFO("Update bsm matching block");
    };

    thdm_parameters::composer.compose_block("WPARAM_MATCH_BSM", src, func);

}