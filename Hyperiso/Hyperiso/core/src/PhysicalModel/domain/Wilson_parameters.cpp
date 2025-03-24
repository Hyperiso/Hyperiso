#include "Wilson_parameters.h"

void WilsonParameters::set_mu(double mu) {
	WilsonParamComposer().update("WPARAM_RUN_SM");

	ParameterProxy wilson_p {ParameterType::WILSON};
	double alphas_mu = wilson_p("WPARAM_RUN_SM", 1);
	double eta_mu = wilson_p("WPARAM_RUN_SM", 2);

	for (int i = 0; i < BWilsonRunningParameters::array_size; ++i) {
        (w_run.etaMuPowers)[i] = std::pow(eta_mu, (BWilsonRunningParameters::ai)[i]);
    }
	for (int i = 0; i < BWilsonRunningParameters::array_size; ++i) {
        (w_run.etaMuPowers2)[i] = std::pow(eta_mu, (BWilsonRunningParameters::ai2)[i]);
    }
	

	LOG_DEBUG("U0,U1, U2 for a scale of " + std::to_string(mu));
	for (int ke = 0; ke < BWilsonRunningParameters::array_size; ++ke) {
        for (int le = 0; le < BWilsonRunningParameters::array_size; ++le) {
            w_run.U0[ke][le] = 0;
            w_run.U1[ke][le] = 0;
            w_run.U2[ke][le] = 0;
			w_run.V0[ke][le] = 0;
			w_run.V1[ke][le] = 0;

            for (int ie = 0; ie < BWilsonRunningParameters::array_size; ++ie) {
                (w_run.U0)[ke][le] += (BWilsonRunningParameters::m00)[ke][le][ie] * (w_run.etaMuPowers)[ie];
                (w_run.U1)[ke][le] += (BWilsonRunningParameters::m10)[ke][le][ie] * (w_run.etaMuPowers)[ie] + (BWilsonRunningParameters::m11)[ke][le][ie] * (w_run.etaMuPowers)[ie] / eta_mu;
                (w_run.U2)[ke][le] += (BWilsonRunningParameters::m20)[ke][le][ie] * (w_run.etaMuPowers)[ie] + (BWilsonRunningParameters::m21)[ke][le][ie] * (w_run.etaMuPowers)[ie] / eta_mu + (BWilsonRunningParameters::m22)[ke][le][ie] * (w_run.etaMuPowers[ie]) / (eta_mu * eta_mu);

				w_run.V0[ke][le]=w_run.V0[ke][le] + BWilsonRunningParameters::l00[ke][le][ie]*pow(eta_mu,BWilsonRunningParameters::ai[ie]);
		        w_run.V1[ke][le]=w_run.V1[ke][le] + BWilsonRunningParameters::l10[ke][le][ie]*pow(eta_mu,BWilsonRunningParameters::ai[ie])+BWilsonRunningParameters::l11[ke][le][ie]*pow(eta_mu,BWilsonRunningParameters::ai[ie]-1.);
            }
			
            LOG_DEBUG("U0[" + std::to_string(ke) + "][" + std::to_string(le) + "]: " + std::to_string((w_run.U0)[ke][le]));
            LOG_DEBUG("U1[" + std::to_string(ke) + "][" + std::to_string(le) + "]: " + std::to_string((w_run.U1)[ke][le]));
            LOG_DEBUG("U2[" + std::to_string(ke) + "][" + std::to_string(le) + "]: " + std::to_string((w_run.U2)[ke][le]));
        }
    }

}

void WilsonParameters::set_mu_W(double mu_W) {
	WilsonParamComposer().update("WPARAM_MATCH_SM");
}
WilsonParameters::WilsonParameters(double mu_W, double mu_h, int gen) {
	this->init_scale_independent_block(gen);
	this->init_matching_block(mu_W);
	this->init_running_block(mu_W, mu_h);
}

void WilsonParameters::init_scale_independent_block(int gen) {
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS"}}};

    auto func = [gen] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        auto xh = std::pow(src.at("MASS")->retrieve(25).get_val() / src.at("MASS")->retrieve(24).get_val(), 2);
		
		int nf = 5;
		int id = 1;
        dep_block->store(id, Parameter({ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id++}}, xh, 0., 0.));
		dep_block->store(id, Parameter({ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id++}}, gen, 0., 0.));
		dep_block->store(id, Parameter({ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id++}}, src.at("MASS")->retrieve(9 + 2 * gen).get_val(), 0., 0.));
		dep_block->store(id, Parameter({ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id++}}, 0.22305, 0., 0.)); //TODO, sw2
		dep_block->store(id, Parameter({ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id++}}, 11.-2./3.*nf, 0., 0.)); //TODO, beta0
    };

    this->composer.compose("WPARAM_SI_SM", src, func);
}

void WilsonParameters::init_matching_block(double mu_W) {
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS", "QCD"}}};

    auto func = [mu_W] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double alphas_muW = QCDHelper::alpha_s(mu_W);
		double mass_top_muW = QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR);
		double mass_b_muW_mbrun = QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR);
		double mass_b_muW_mbpole = QCDHelper::msbar_mass(5, mu_W, MassType::POLE);
		double mass_c_muW = QCDHelper::msbar_mass(4, mu_W, MassType::POLE);;
		double xt = std::pow(mass_top_muW / src.at("MASS")->retrieve(24).get_val(), 2);

		dep_block->store(1, Parameter({ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 1}}, alphas_muW, 0., 0.));
		dep_block->store(LhaID(2, 1), Parameter({ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 1)}}, xt, 0., 0.));
		dep_block->store(LhaID(2, 2), Parameter({ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 2)}}, std::pow(xt, 2), 0., 0.));
		dep_block->store(LhaID(2, 3), Parameter({ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 3)}}, std::pow(xt, 3), 0., 0.));
		dep_block->store(LhaID(2, 4), Parameter({ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(2, 4)}}, std::pow(xt, 4), 0., 0.));
		dep_block->store(4, Parameter({ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 4}}, mass_c_muW, 0., 0.));
		dep_block->store(LhaID(5, 1), Parameter({ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}}, mass_b_muW_mbrun, 0., 0.));
		dep_block->store(LhaID(5, 2), Parameter({ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 2)}}, mass_b_muW_mbpole, 0., 0.));
		dep_block->store(6, Parameter({ParamId{ParameterType::WILSON, "WPARAM_MATCH_SM", 6}}, mass_top_muW, 0., 0.));
    };

    this->composer.compose("WPARAM_MATCH_SM", src, func);
}

void WilsonParameters::init_running_block(double mu_W, double mu_h) {
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS", "QCD"}}};

    auto func = [mu_W, mu_h] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double alphas_mu = QCDHelper::alpha_s(mu_h);	
		double eta = QCDHelper::alpha_s(mu_W) / alphas_mu;

		std::array<double, BWilsonRunningParameters::array_size> etaMuPowers = {};
		std::array<double, BWilsonRunningParameters::array_size> etaMuPowers2 = {};
		for (int i = 0; i < BWilsonRunningParameters::array_size; ++i) {
			(etaMuPowers)[i] = std::pow(eta, (BWilsonRunningParameters::ai)[i]);
		}
		for (int i = 0; i < BWilsonRunningParameters::array_size; ++i) {
			(etaMuPowers2)[i] = std::pow(eta, (BWilsonRunningParameters::ai2)[i]);
		}

		std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> U0 = {};
		std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> U1 = {};
		std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> U2 = {};
		std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> V0 = {};
		std::array<std::array<double, BWilsonRunningParameters::array_size>, BWilsonRunningParameters::array_size> V1 = {};
		
		dep_block->store(1, Parameter({ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 1}}, alphas_mu, 0., 0.));
		dep_block->store(2, Parameter({ParamId{ParameterType::WILSON, "WPARAM_RUN_SM", 2}}, eta, 0., 0.));
    };

    this->composer.compose("WPARAM_RUN_SM", src, func);
}
