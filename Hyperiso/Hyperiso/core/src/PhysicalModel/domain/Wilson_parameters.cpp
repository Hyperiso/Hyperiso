#include "Wilson_parameters.h"

// void Wilson_parameters::SetMu(double mu) {
// 	alphas_mu=QCDHelper::alpha_s(mu);	
// 	eta_mu=alphas_muW/alphas_mu;

// 	for (int i = 0; i < array_size; ++i) {
//         (etaMuPowers)[i] = std::pow(eta_mu, (ai)[i]);
//     }
// 	for (int i = 0; i < array_size; ++i) {
//         (etaMuPowers2)[i] = std::pow(eta_mu, (ai2)[i]);
//     }
	

// 	LOG_DEBUG("U0,U1, U2 for a scale of " + std::to_string(mu));
// 	for (int ke = 0; ke < array_size; ++ke) {
//         for (int le = 0; le < array_size; ++le) {
//             (U0)[ke][le] =0;
//             (U1)[ke][le] = 0;
//             (U2)[ke][le] = 0;
// 			V0[ke][le] = 0;
// 			V1[ke][le] = 0;
//             for (int ie = 0; ie < array_size; ++ie) {
//                 (U0)[ke][le] += (m00)[ke][le][ie] * (etaMuPowers)[ie];
//                 (U1)[ke][le] += (m10)[ke][le][ie] * (etaMuPowers)[ie] + (m11)[ke][le][ie] * (etaMuPowers)[ie] / eta_mu;
//                 (U2)[ke][le] += (m20)[ke][le][ie] * (etaMuPowers)[ie] + (m21)[ke][le][ie] * (etaMuPowers)[ie] / eta_mu + (m22)[ke][le][ie] * (etaMuPowers[ie]) / (eta_mu * eta_mu);

// 				V0[ke][le]=V0[ke][le] + l00[ke][le][ie]*pow(eta_mu,ai[ie]);
// 		        V1[ke][le]=V1[ke][le] + l10[ke][le][ie]*pow(eta_mu,ai[ie])+l11[ke][le][ie]*pow(eta_mu,ai[ie]-1.);
//             }
			
//             LOG_DEBUG("U0[" + std::to_string(ke) + "][" + std::to_string(le) + "]: " + std::to_string((U0)[ke][le]));
//             LOG_DEBUG("U1[" + std::to_string(ke) + "][" + std::to_string(le) + "]: " + std::to_string((U1)[ke][le]));
//             LOG_DEBUG("U2[" + std::to_string(ke) + "][" + std::to_string(le) + "]: " + std::to_string((U2)[ke][le]));
//         }
//     }

// }

WilsonParameters::WilsonParameters(double mu_W, double mu_h, int gen) {
	this->init_scale_independent_block(gen);
	this->init_matching_block(mu_W);
	this->init_running_block(mu_W, mu_h);
}

void WilsonParameters::init_scale_independent_block(int gen) {
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS"}}};

    auto func = [gen] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        auto xh = std::pow(src.at("MASS")->retrieve(25).get_val() / src.at("MASS")->retrieve(24).get_val(), 2);

		int id = 1;
        dep_block->store(id, Parameter({ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id++}}, xh, 0., 0.));
		dep_block->store(id, Parameter({ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id++}}, gen, 0., 0.));
		dep_block->store(id, Parameter({ParamId{ParameterType::WILSON, "WPARAM_SI_SM", id++}}, src.at("MASS")->retrieve(9 + 2 * gen).get_val(), 0., 0.));
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
		
    };

    this->composer.compose("WPARAM_RUN_SM", src, func);
}
