#include "Wilson_parameters.h"

void WilsonParameterHelper::init(int gen) {
	if (WilsonParameterHelper::initialized) {
		return;
	}

	LOG_DEBUG("Initializing WilsonParameterHelper");
	WilsonParameterHelper::init_scale_independent_block(gen);
	WilsonParameterHelper::init_matching_block();
	WilsonParameterHelper::init_running_block();
	WilsonParameterHelper::initialized = true;
}

void WilsonParameterHelper::init_scale_independent_block(int gen) {
	LOG_DEBUG("Init scale-independent wparam block");
	ParameterProxy(ParameterType::SM);
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"SMINPUTS", "MASS"}}};

    auto func = [gen] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double xh = pow(src.at("MASS")->retrieve(25)->get_val() / src.at("MASS")->retrieve(24)->get_val(), 2);
		
		int nf = 5;
		int id = 1;
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, xh, 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, gen, 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, src.at("MASS")->retrieve(9 + 2 * gen)->get_val(), 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, src.at("SMINPUTS")->retrieve({7, 1})->get_val(), 0., 0.));
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id}, 11.-2./3.*nf, 0., 0.)); //TODO, beta0
    };
    WilsonParameterHelper::composer.compose_block("WPARAM_SI_SM", src, func);
}

void WilsonParameterHelper::init_matching_block() {
	LOG_DEBUG("Init matching scale dependent wparam block");
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS", "QCD"}}, {ParameterType::WILSON, {"EW_SCALE"}}};

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double mu_W = src.at("EW_SCALE")->retrieve(1)->get_val();
        double alphas_muW = QCDHelper::alpha_s(mu_W);
		double mass_top_muW = QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR);
		double mass_b_muW_mbrun = QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR);
		double mass_b_muW_mbpole = QCDHelper::msbar_mass(5, mu_W, MassType::POLE);
		double mass_c_muW = QCDHelper::msbar_mass(4, mu_W, MassType::POLE);

		double m_W = src.at("MASS")->retrieve(24)->get_val();
		double xt = pow(mass_top_muW / m_W, 2);
		double L = log(std::pow(mu_W / m_W, 2));
		double xtW = pow(QCDHelper::msbar_mass(6, m_W) / m_W, 2); // mass top at pole for mtot param
		double xtt = pow(src.at("QCD")->retrieve(6)->get_val() / m_W, 2.); // 24 -> W

		dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 1}, alphas_muW, 0., 0.));
		dep_block->store_or_assign(LhaID(2, 1), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)}, xt, 0., 0.));
		dep_block->store_or_assign(LhaID(2, 2), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 2)}, std::pow(xt, 2), 0., 0.));
		dep_block->store_or_assign(LhaID(2, 3), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 3)}, std::pow(xt, 3), 0., 0.));
		dep_block->store_or_assign(LhaID(2, 4), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 4)}, std::pow(xt, 4), 0., 0.));
		dep_block->store_or_assign(3, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 3}, L, 0., 0.));
		dep_block->store_or_assign(4, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 4}, mass_c_muW, 0., 0.));
		dep_block->store_or_assign(LhaID(5, 1), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, mass_b_muW_mbrun, 0., 0.));
		dep_block->store_or_assign(LhaID(5, 2), std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 2)}, mass_b_muW_mbpole, 0., 0.));
		dep_block->store_or_assign(6, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 6}, mass_top_muW, 0., 0.));
		dep_block->store_or_assign(7, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 7}, xtW, 0., 0.));
		dep_block->store_or_assign(8, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 8}, xtt, 0., 0.));

		LOG_DEBUG("Update matching block");
    };

    WilsonParameterHelper::composer.compose_block("WPARAM_MATCH_SM", src, func);
}

void WilsonParameterHelper::init_running_block() {
	LOG_DEBUG("Init running scale dependent wparam block");
	std::unordered_map<ParameterType, std::vector<std::string>> src = {
		{ParameterType::WILSON, {"EW_SCALE", "B_SCALE"}}
	};

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
		double mu_W = src.at("EW_SCALE")->retrieve(1)->get_val();
		double mu_h = src.at("B_SCALE")->retrieve(1)->get_val();
        double alphas_mu = QCDHelper::alpha_s(mu_h);	
		double eta = QCDHelper::alpha_s(mu_W) / alphas_mu;
		
		dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 1}, alphas_mu, 0., 0.));
		dep_block->store_or_assign(2, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 2}, eta, 0., 0.));
    };

    WilsonParameterHelper::composer.compose_block("WPARAM_RUN_SM", src, func);
}

void WilsonParameterHelper::cleanup() {
	WilsonParameterHelper::initialized = false;
}