#include "Wilson_parameters.h"

void WilsonParameterHelper::set_mu(double mu) {
	if (fpeq(mu, WilsonParameterHelper::current_mu_h)) {
		return;
	}

	WilsonParamComposer().update("WPARAM_RUN_SM");
	WilsonParameterHelper::current_mu_h = mu;
}

void WilsonParameterHelper::set_mu_W(double mu_W) {
	if (fpeq(mu_W, WilsonParameterHelper::current_mu_W)) {
		std::cout << "WilsonParametersHelper already set a scale " << mu_W << std::endl;
		return;
	}
	std::cout << "mmh" << std::endl;
    WilsonParamComposer().update("WPARAM_MATCH_SM");
	WilsonParamComposer().update("WPARAM_RUN_SM");
	WilsonParameterHelper::current_mu_W = mu_W;
}

void WilsonParameterHelper::init(double mu_W, double mu_h, int gen) {
	if (WilsonParameterHelper::initialized) {
		return;
	}

	std::cout << "Initializing WilsonParameterHelper at scales " << mu_W << " and " << mu_h << std::endl;
	WilsonParameterHelper::current_mu_W = mu_W;
	WilsonParameterHelper::current_mu_h = mu_h;
	WilsonParameterHelper::init_scale_independent_block(gen);
	WilsonParameterHelper::init_matching_block(mu_W);
	WilsonParameterHelper::init_running_block(mu_W, mu_h);
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

void WilsonParameterHelper::init_matching_block(double mu_W) {
	LOG_DEBUG("Init matching scale dependent wparam block at mu_W =", mu_W);
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS" /*, "QCD"*/}}};

    auto func = [mu_W] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double alphas_muW = QCDHelper::alpha_s(mu_W);
		double mass_top_muW = QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR);
		double mass_b_muW_mbrun = QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR);
		double mass_b_muW_mbpole = QCDHelper::msbar_mass(5, mu_W, MassType::POLE);
		double mass_c_muW = QCDHelper::msbar_mass(4, mu_W, MassType::POLE);

		double m_W = src.at("MASS")->retrieve(24)->get_val();
		double xt = pow(mass_top_muW / m_W, 2);
		double L = log(std::pow(mu_W / m_W, 2));
		double xtW = pow(QCDHelper::msbar_mass(6, m_W) / m_W, 2); // mass top at pole for mtot param
		double xtt = pow(QCDHelper::mass_t_msbar() / m_W, 2.); // 24 -> W

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

		LOG_INFO("Update matching block");
    };

    WilsonParameterHelper::composer.compose_block("WPARAM_MATCH_SM", src, func);
}

void WilsonParameterHelper::init_running_block(double mu_W, double mu_h) {
	LOG_DEBUG("Init running scale dependent wparam block at mu_W =", mu_W, "and mu_h = ", mu_h);
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS"/*, "QCD" */}}};

    auto func = [mu_W, mu_h] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double alphas_mu = QCDHelper::alpha_s(mu_h);	
		double eta = QCDHelper::alpha_s(mu_W) / alphas_mu;
		
		dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 1}, alphas_mu, 0., 0.));
		dep_block->store_or_assign(2, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 2}, eta, 0., 0.));
    };

    WilsonParameterHelper::composer.compose_block("WPARAM_RUN_SM", src, func);
}