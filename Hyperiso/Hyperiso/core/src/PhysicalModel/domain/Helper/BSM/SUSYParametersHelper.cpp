#include "SUSYParametersHelper.h"

void susy_parameters::init(int gen, WGroupId grp) {

	if (initialized) {
		return;
	}

	LOG_INFO("Initializing scale independent SUSY parameters");
	init_scale_independent_block(gen);
	LOG_INFO("Initializing matching SUSY parameters");
	init_matching_block();
	LOG_INFO("Initializing epsilon block");
	init_epsilon_block();
	LOG_INFO("Done");

	initialized = true;
}

void susy_parameters::init_epsilon_block() {
    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::SM, {"MASS", "SMINPUTS"}}, 
        {ParameterType::BSM, {"MASS", "GAUGE", "HMIX", "MSOFT", "AD", "AU", "YD", "YU", "SBOTMIX", "STOPMIX", "UMIX", "VMIX", "NMIX", "ALPHA"}},
        {ParameterType::WILSON, {"WPARAM_SI_SM"}}
    };

    auto func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {

        src.get_val("ALPHA", {});

        double g2 = src.get_val("GAUGE", 2);
        double alpha_em = src.get_val("SMINPUTS", 1);

        double m_ds = src.get_val("MASS", 1000001);
        double m_us = src.get_val("MASS", 1000002);
        double m_ss = src.get_val("MASS", 1000003);
        double m_cs = src.get_val("MASS", 1000004);
        double m_bs = src.get_val("MASS", 1000005);
        double m_ts = src.get_val("MASS", 1000006);

        double m_d2s = src.get_val("MASS", 2000001);
        double m_u2s = src.get_val("MASS", 2000002);
        double m_s2s = src.get_val("MASS", 2000003);
        double m_c2s = src.get_val("MASS", 2000004);
        double m_b2s = src.get_val("MASS", 2000005);
        double m_t2s = src.get_val("MASS", 2000006);

        double m_gluino = src.get_val("MASS", 1000021);
        double m_c1 = src.get_val("MASS", 1000024);
        double m_c2 = src.get_val("MASS", 1000037);
        double mu_Q = src.get_val("HMIX", 1);

        double mqL3 = src.get_val("MSOFT", 43);
        double mbR = src.get_val("MSOFT", 49);

        double m_G = src.get_val("MASS", 1000039);

        std::map<int,int> neutralino = {{0, 1000022},{1, 1000023},{2, 1000025},{3, 1000035}};


        std::vector<double> m_neutralino = {src.get_val("MASS", neutralino[0]), src.get_val("MASS", neutralino[1]), src.get_val("MASS", neutralino[2]), src.get_val("MASS", neutralino[3])};
        double ad_22 = src.get_val("AD", {3, 3});
        double au_22 = src.get_val("AU", {3, 3});

        double yu_22 = src.get_val("YU", {3, 3});
        double yd_22 = src.get_val("YD", {3, 3});

        double sbot_mix_00 = src.get_val("SBOTMIX", {0+1, 0+1});
        double sbot_mix_01 = src.get_val("SBOTMIX", {0+1, 1+1});

        double stop_mix_00 = src.get_val("STOPMIX", {0+1, 0+1});
        double stop_mix_01 = src.get_val("STOPMIX", {0+1, 1+1});

        double umix_01 = src.get_val("UMIX", {0+1, 1+1});
        double umix_11 = src.get_val("UMIX", {1+1, 1+1});

        double vmix_01 = src.get_val("VMIX", {0+1, 1+1});
        double vmix_11 = src.get_val("VMIX", {1+1, 1+1});
        
        double sw2 = src.get_val("WPARAM_SI_SM", 4);

        double alphas_MSOFT = QCDHelper::alpha_s(2.448e3); //TODO better
        double MSOFT = ParameterProxy(ParameterType::BSM).get_scale("MSOFT"); //TODO : better things to do with scale

        std::cout << "MSOFT" << MSOFT << std::endl; 
        double tan_beta = src.get_val("HMIX", 2);

        double factor = 2.0 / 3.0 * alphas_MSOFT / M_PI;

        double M_2 = src.get_val("MSOFT", 2);


        double term1 =  (ad_22 / tan_beta - mu_Q) / m_gluino *
                H2(m_bs * m_bs / m_gluino / m_gluino, m_b2s * m_b2s / m_gluino / m_gluino);
        double term2 = -0.5 * (B(m_gluino, m_bs, MSOFT) + B(m_gluino, m_b2s, MSOFT)) / tan_beta;
        double term3 = 1.0 / alpha_em / sw2 / 4.0 / M_PI * (mu_Q * M_2) * 
                (sbot_mix_00 * sbot_mix_00 * H2(M_2 * M_2 / m_bs / m_bs, mu_Q * mu_Q / m_bs / m_bs) / m_bs / m_bs / 2.0 +
                sbot_mix_01 * sbot_mix_01 * H2(M_2 * M_2 / m_b2s / m_b2s, mu_Q * mu_Q / m_b2s / m_b2s) / m_b2s / m_b2s / 2.0);


        double epsilon_0 = factor * (term1 + term2) + term3;


        term1 = yu_22 * yu_22 / 16.0 / M_PI / M_PI * 
                    (mu_Q / tan_beta - au_22) * 
                    ((umix_01 * vmix_01 / m_c1 * 
                        H2(m_ts * m_ts / m_c1 / m_c1, m_t2s * m_t2s / m_c1 / m_c1)) +
                        (umix_11 * vmix_11 / m_c2 * 
                        H2(m_ts * m_ts / m_c2 / m_c2, m_t2s * m_t2s / m_c2 / m_c2)));
        
        term2 = 1.0 / alpha_em / sw2 / 4.0 / M_PI * (mu_Q * M_2) * 
                    ((stop_mix_00 * stop_mix_00 * 
                        H2(M_2 * M_2 / m_ts / m_ts, mu_Q * mu_Q / m_ts / m_ts) / m_ts / m_ts) +
                        (stop_mix_01 * stop_mix_01* 
                        H2(M_2 * M_2 / m_t2s / m_t2s, mu_Q * mu_Q / m_t2s / m_t2s) / m_t2s / m_t2s));

        double epsilon_2 = term1 + term2;


        //b
        double epsilon_b = epsilon_0 + epsilon_2;

        //bp
        int nb_neut = (m_G == 0.) ? 4 : 5; //mass_neut[5] is gravitino ?

        
        double epsilonbp = 2.0 / 3.0 * alphas_MSOFT / M_PI * 
                        (ad_22/ tan_beta - mu_Q) / m_gluino * 
                        (stop_mix_00 * stop_mix_00 * sbot_mix_00 * sbot_mix_00*
                            H2(m_ts * m_ts / m_gluino / m_gluino, m_b2s * m_b2s / m_gluino / m_gluino) +
                            stop_mix_00 * stop_mix_00 * sbot_mix_01 * sbot_mix_01 *
                            H2(m_ts * m_ts / m_gluino / m_gluino, m_bs * m_bs / m_gluino / m_gluino) +
                            stop_mix_01 * stop_mix_01 * sbot_mix_00 * sbot_mix_00 *
                            H2(m_t2s * m_t2s / m_gluino / m_gluino, m_b2s * m_b2s / m_gluino / m_gluino) +
                            stop_mix_01 * stop_mix_01 * sbot_mix_01 * sbot_mix_01 *
                            H2(m_t2s * m_t2s / m_gluino / m_gluino, m_bs * m_bs / m_gluino / m_gluino));

        for(int ie = 0; ie < nb_neut; ++ie) {
            epsilonbp += yu_22 * yu_22 / 16.0 / M_PI / M_PI * 
                        src.get_val("NMIX", {ie + 1, 3 + 1}) * src.get_val("NMIX", {ie + 1, 2 + 1}) * 
                        (au_22 - mu_Q / tan_beta) / m_neutralino[ie] *
                        (stop_mix_00 * stop_mix_00 * sbot_mix_00 * sbot_mix_00 *
                        H2(m_t2s * m_t2s / m_neutralino[ie] / m_neutralino[ie], m_bs * m_bs / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_00 * stop_mix_00 * sbot_mix_01 * sbot_mix_01 *
                        H2(m_t2s * m_t2s / m_neutralino[ie] / m_neutralino[ie], m_b2s * m_b2s / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_01 * stop_mix_01 * sbot_mix_00 * sbot_mix_00 *
                        H2(m_ts * m_ts / m_neutralino[ie] / m_neutralino[ie], m_bs * m_bs / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_01* stop_mix_01 * sbot_mix_01 * sbot_mix_01 *
                        H2(m_ts * m_ts / m_neutralino[ie] / m_neutralino[ie], m_b2s * m_b2s / m_neutralino[ie] / m_neutralino[ie]));
        }

        epsilonbp += 1.0 / alpha_em / sw2 / 4.0 / M_PI * 
                    (mu_Q * M_2) * 
                    ((stop_mix_00 * stop_mix_00 *
                    H2(M_2 * M_2 / m_ts / m_ts, mu_Q * mu_Q / m_ts / m_ts) / m_ts / m_ts +
                    stop_mix_01 * stop_mix_01 *
                    H2(M_2 * M_2 / m_t2s / m_t2s, mu_Q * mu_Q / m_t2s / m_t2s) / m_t2s / m_t2s) / 2.0 +
                    (sbot_mix_00 * sbot_mix_00 *
                    H2(M_2 * M_2 / m_bs / m_bs, mu_Q * mu_Q / m_bs / m_bs) / m_bs / m_bs +
                    sbot_mix_01 * sbot_mix_01 *
                    H2(M_2 * M_2 / m_b2s / m_b2s, mu_Q * mu_Q / m_b2s / m_b2s) / m_b2s / m_b2s));


        double epsilon0p = -2.0 / 3.0 * alphas_MSOFT / M_PI * 
                        (mu_Q + au_22 / tan_beta) / m_gluino *
                        (stop_mix_00 * stop_mix_00 * 
                            H2(m_t2s * m_t2s / m_gluino / m_gluino, m_ss * m_ss / m_gluino / m_gluino) +
                            stop_mix_01 * stop_mix_01 * 
                            H2(m_ts * m_ts / m_gluino / m_gluino, m_ss * m_ss / m_gluino / m_gluino));

        for(int ie = 0; ie < nb_neut; ++ie) {
            epsilon0p += yd_22 * yd_22 / 16.0 / M_PI / M_PI * 
                        src.get_val("NMIX", {ie+1, 3+1}) * src.get_val("NMIX", {ie+1, 2+1}) * 
                        (mu_Q / tan_beta) / m_neutralino[ie] *
                        (stop_mix_00 * stop_mix_00 * sbot_mix_00 * sbot_mix_00 * 
                        H2(m_ts * m_ts / m_neutralino[ie] / m_neutralino[ie], m_b2s * m_b2s / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_00 * stop_mix_00 * sbot_mix_01 * sbot_mix_01 * 
                        H2(m_ts * m_ts / m_neutralino[ie] / m_neutralino[ie], m_bs * m_bs / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_01 * stop_mix_01 * sbot_mix_00 * sbot_mix_00 * 
                        H2(m_t2s * m_t2s / m_neutralino[ie] / m_neutralino[ie], m_b2s * m_b2s / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_01 * stop_mix_01 * sbot_mix_01* sbot_mix_01 * 
                        H2(m_t2s * m_t2s / m_neutralino[ie] / m_neutralino[ie], m_bs * m_bs / m_neutralino[ie] / m_neutralino[ie]));

        }

        term1 = 1.0 / 16.0 / M_PI / M_PI * 
        (yd_22 * yd_22 * ad_22 / mu_Q * 
        H2(std::pow(mqL3 / mu_Q, 2), std::pow(mbR / mu_Q, 2))); 

        term2 = -g2 * g2 * M_2 / mu_Q * 
            H2(std::pow(mqL3 / mu_Q, 2), std::pow(M_2 / mu_Q, 2)) / 16.0 / M_PI / M_PI;

        double epsilon_1p = term1 + term2;
        
        double epsfac=pow((1.+epsilon_b*tan_beta),2.);

        dep_block->store_or_assign({0,1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", {0,1}}, epsilon_0, 0., 0.));
		dep_block->store_or_assign({0,2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", {0,2}}, epsilon0p, 0., 0.));
		dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 1}, epsilon_1p, 0., 0.));
		dep_block->store_or_assign(2, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 2}, epsilon_2, 0., 0.));
		dep_block->store_or_assign(3, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 3}, epsilon_b, 0., 0.));
        dep_block->store_or_assign(4, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 4}, epsilonbp, 0., 0.));
        dep_block->store_or_assign(5, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 5}, epsfac, 0., 0.)); 

    };

    iblock_c->compose_block("EPSILON_SUSY", src, func);
}

void susy_parameters::init_scale_independent_block(int gen) {

	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS", "GAUGE", "VCKM"}}, {ParameterType::BSM, {"MASS", "HMIX", "STOPMIX"}}};

    auto func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
        double mW = src.get_val("MASS", 24);
		double alphas_mg = QCDHelper::alpha_s(src.get_val("MASS", 1000021));
		double ag = 1.0 - 7.0 / (12.0 * Pi) * alphas_mg;
		double aY = 1.0 + alphas_mg / (4.0 * Pi);

		double kappa = 1.0 / (pow(src.get_val("GAUGE", 2), 2.) * 
						std::real((src.get_val("VCKM", {2,2}))*(src.get_val("VCKM", {2,1})))); //VCKM 33 et 32
		//TODO : careful with real

		double kappaFactor = -0.5 * kappa;

		double tanb = src.get_val("HMIX", 2);
		double beta = std::atan(tanb);
		double z = pow(src.get_val("MASS", 37) / mW, 2.);
		double sinb = std::sin(std::atan(tanb));
		double cosb = std::cos(std::atan(tanb));
		double ct = src.get_val("STOPMIX", {2,2});
		double st = src.get_val("STOPMIX", {1,2});

        double lu = 1./tanb;
        double ld = -tanb;

		Array1D_4 ME = {src.get_val("MASS", 11), src.get_val("MASS", 13), src.get_val("MASS", 15)}; 

		Array1D_3 Mch = {src.get_val("MASS", 1000024), src.get_val("MASS", 1000037)};

		Array1D_7 MsqU = {src.get_val("MASS", 1000002), src.get_val("MASS", 1000004), src.get_val("MASS", 1000006), 
				src.get_val("MASS", 2000002), src.get_val("MASS", 2000004), src.get_val("MASS", 2000006)};

		Array1D_7 MsqD = {src.get_val("MASS", 1000001), src.get_val("MASS", 1000003), src.get_val("MASS", 1000005), 
				src.get_val("MASS", 2000001), src.get_val("MASS", 2000003), src.get_val("MASS", 2000005)};

		Array1D_4 Msn = {src.get_val("MASS", 1000012), src.get_val("MASS", 1000014), src.get_val("MASS", 1000016)};
		
		const size_t NumSquarks = 6;
		Array2D_7x7 sU_mix; //ERROR
		bool isNonZeroMix = true;
		for (size_t i = 0; i < NumSquarks; ++i) {
			double product = 1.0;
			for (size_t j = 0; j < 6; ++j) {
				product *= sU_mix[i][j]; //TODO, ISSUE
			}
			if (product == 0.0) {
				isNonZeroMix = false;
				break;
			}
		}

		if (isNonZeroMix) {
			printf("eheheheh\n");
			std::sort(MsqU.begin(), MsqU.end());
		}
        int id {1};
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, z, 0., 0.)); //1
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, cosb, 0., 0.)); //2
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, sinb, 0., 0.)); //3
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, ct, 0., 0.)); //4
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, st, 0., 0.)); //5
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, beta, 0., 0.)); //6
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, lu, 0., 0.)); //7
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, ld, 0., 0.)); //8
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, alphas_mg, 0., 0.)); //9
        dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, ag, 0., 0.)); //10
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, aY, 0., 0.)); //11
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, ME[0], 0., 0.)); //12
		dep_block->store_or_assign({id, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 1}}, ME[1], 0., 0.)); //12
		dep_block->store_or_assign({id++, 2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 2}}, ME[2], 0., 0.)); //12
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, Mch[0], 0., 0.)); //13
		dep_block->store_or_assign({id++, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 1}}, Mch[1], 0., 0.)); //13
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqU[0], 0., 0.)); //14
		dep_block->store_or_assign({id, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 1}}, MsqU[1], 0., 0.)); //14
		dep_block->store_or_assign({id, 2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 2}}, MsqU[2], 0., 0.)); //14
		dep_block->store_or_assign({id, 3}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 3}}, MsqU[3], 0., 0.)); //14
		dep_block->store_or_assign({id, 4}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 4}}, MsqU[4], 0., 0.)); //14
		dep_block->store_or_assign({id++, 5}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 5}}, MsqU[5], 0., 0.)); //14
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, MsqD[0], 0., 0.)); //15
		dep_block->store_or_assign({id, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 1}}, MsqD[1], 0., 0.)); //15
		dep_block->store_or_assign({id, 2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 2}}, MsqD[2], 0., 0.)); //15
		dep_block->store_or_assign({id, 3}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 3}}, MsqD[3], 0., 0.)); //15
		dep_block->store_or_assign({id, 4}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 4}}, MsqD[4], 0., 0.)); //15
		dep_block->store_or_assign({id++, 5}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 5}}, MsqD[5], 0., 0.)); //15
		dep_block->store_or_assign({id, 0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM",{id, 0}}, Msn[0], 0., 0.)); //16
		dep_block->store_or_assign({id, 1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 1}}, Msn[1], 0., 0.)); //16
		dep_block->store_or_assign({id++, 2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", {id, 2}}, Msn[2], 0., 0.)); //16
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, (double)isNonZeroMix, 0., 0.)); //17
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, kappaFactor, 0., 0.)); //18
		dep_block->store_or_assign(id++, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_SI_BSM", id}, kappa, 0., 0.)); //19
    };

    iblock_c->compose_block("WPARAM_SI_BSM", src, func);


}

void susy_parameters::init_matching_block() {
	std::unordered_map<ParameterType, std::vector<std::string>> src = {{ParameterType::SM, {"MASS"}}, {ParameterType::BSM, {"MASS"}},
																	{ParameterType::WILSON, {"WPARAM_MATCH_SM"}}};

    auto func = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {
        double yt= pow(src.get_val("WPARAM_MATCH_SM", 6)/src.get_val("MASS", 37),2.); // param->mass_H (25)
		Array1D_4 MU = {src.get_val("MASS", 2), src.get_val("MASS", 4), src.get_val("WPARAM_MATCH_SM", 6)}; //TODO : size 3 not 4
		Array1D_4 MD = {src.get_val("MASS", 2), src.get_val("MASS", 3), src.get_val("WPARAM_MATCH_SM", {5,1})}; //TODO : size 3 not 4 // TODO : MD[0] -> mu like superiso but why ?

        dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", 1}, yt, 0., 0.));
		dep_block->store_or_assign({2,0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,0}}, MU[0], 0., 0.));
		dep_block->store_or_assign({2,1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,1}}, MU[1], 0., 0.));
		dep_block->store_or_assign({2,2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {2,2}}, MU[2], 0., 0.));
		dep_block->store_or_assign({3,0}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {3,0}}, MD[0], 0., 0.));
		dep_block->store_or_assign({3,1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {3,1}}, MD[1], 0., 0.));
		dep_block->store_or_assign({3,2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "WPARAM_MATCH_BSM", {3,2}}, MD[2], 0., 0.));

    };

    iblock_c->compose_block("WPARAM_MATCH_BSM", src, func);

	std::unordered_map<ParameterType, std::vector<std::string>> src_matrix = {{ParameterType::SM, {"MASS", "VCKM", "GAUGE"}}, {ParameterType::BSM, {"UMIX", "VMIX"}},
																				{ParameterType::WILSON, {"WPARAM_SI_BSM", "WPARAM_MATCH_SM", "WPARAM_MATCH_BSM"}}};

    auto func_matrix = [] (const BlockSrc& src, std::shared_ptr<DependentBlock> dep_block) {

		Array2D_7x4 Gamma_UL {};
		Array2D_7x4 Gamma_UR {};
		Array2D_4x4 Gamma_NL {};
		Array2D_4x4 Gamma_NR {};
		Array2D_7x7 Gamma_U{};
		Array2D_7x7 I_LR{};
		Array2D_7x7 P_U{};
		Array3D_3x7x4 X_UL{};
		Array3D_3x7x4 X_UR{};
		Array3D_3x7x4 X_NL{};
		Array3D_3x7x4 X_NR{};
		std::array<std::array<std::array<std::array<double, 4>, 4>, 3>, 7> G_aimn;

		complex_t c11 = src.get_val("VCKM", {0,0});
		complex_t c12 = src.get_val("VCKM", {0,1});
		complex_t c13 = src.get_val("VCKM", {0,2});
		complex_t c21 = src.get_val("VCKM", {1,0});
		complex_t c22 = src.get_val("VCKM", {1,1});
		complex_t c23 = src.get_val("VCKM", {1,2});
		complex_t c31 = src.get_val("VCKM", {2,0});
		complex_t c32 = src.get_val("VCKM", {2,1});
		complex_t c33 = src.get_val("VCKM", {2,2});

		complex_t complexTerm = -(c32 * c33 + c22 * c23) / c12;

		Array2D_4x4_I VCKM = {{
			{
				c11,
				c12,
				c13
			},
			{
				c21,
				c22,
				c23
			},
			{
				c31,
				c32,
				c33
			}
		}};

		double B0c1 = 0.0, B0c2 = 0.0, B90c = 0.0, B100c = 0.0, C90c = 0.0, D90c = 0.0;
    	bool test;

        double mW = src.get_val("MASS", 24);
		double g2 = src.get_val("GAUGE", 2);
		if (src.get_val("WPARAM_SI_BSM", 17)) {
			std::cout << "SHOULD NOT BE HERE" << std::endl;
			Array2D_7x7 sU_mix; //TODO : wtf
			const size_t NumSquarks = 6;
			for (size_t ae = 0; ae < NumSquarks; ++ae) {
				for (size_t ie = 0; ie < 3; ++ie) {
					Gamma_UL[ae][ie] = sU_mix[ae][ie];
					Gamma_UR[ae][ie] = sU_mix[ae][ie + 3];
				}
			}
	}
		else {
			Gamma_UL[0][0] = 1.0; 
			Gamma_UL[1][1] = 1.0;
			Gamma_UL[2][2] = src.get_val("WPARAM_SI_BSM", 4);
			Gamma_UL[5][2] = -src.get_val("WPARAM_SI_BSM", 5);

			Gamma_UR[3][0] = 1.0;
			Gamma_UR[4][1] = 1.0;
			Gamma_UR[2][2] = src.get_val("WPARAM_SI_BSM", 5);
			Gamma_UR[5][2] = src.get_val("WPARAM_SI_BSM", 4);
		}

		for (int ae = 0; ae < 6; ++ae) {
			for (int ie = 0; ie < 3; ++ie) {
				Gamma_U[ae][ie] = Gamma_UL[ae][ie];
				Gamma_U[ae][ie+3] = Gamma_UR[ae][ie];
				if (ae <3 && ae==ie) {
					Gamma_NL[ae][ie] = 1.;
					Gamma_NR[ae][ie] = 1.;
				}
			}
		}

		I_LR.fill({});
		for (int i = 0; i < 3; ++i) {
			I_LR[i][i] = 1.;
			I_LR[i+3][i+3] = -1.;
		}

		for (int ae = 0; ae < 6; ++ae) {
			for (int be = 0; be < 6; ++be) {
				for (int ce = 0; ce < 6; ++ce) {
					for (int de = 0; de < 6; ++de) {
						P_U[ae][be] = Gamma_U[ae][ce] * I_LR[ce][de] * Gamma_U[be][de];
					}
				}
			}
			
		}
		
		for (int ie = 0; ie < 2; ++ie) {
			for (int ae = 0; ae < 6; ++ae) {
				for (int be = 0; be < 3; ++be) {
					X_UL[ie][ae][be] = 0.0;
					X_UR[ie][ae][be] = 0.0;

					for (int ce = 0; ce < 3; ++ce) {
						X_UL[ie][ae][be] += -g2 * (
							src.get_val("WPARAM_SI_BSM", 10) * src.get_val("VMIX", {ie+1, 0+1}) * Gamma_UL[ae][ce] -
							src.get_val("WPARAM_SI_BSM", 11) * src.get_val("VMIX", {ie+1, 1+1}) * Gamma_UR[ae][ce] * src.get_val("WPARAM_MATCH_BSM", {2,ce}) / (sqrt(2.0) * mW * src.get_val("WPARAM_SI_BSM", 3))
						) * std::real(VCKM[ce][be]);
						X_UR[ie][ae][be] += g2 * src.get_val("WPARAM_SI_BSM", 11) * src.get_val("UMIX", {ie+1, 1+1}) * Gamma_UL[ae][ce] * std::real(VCKM[ce][be]) * src.get_val("WPARAM_MATCH_BSM", {3,be}) / (sqrt(2.0) * mW * src.get_val("WPARAM_SI_BSM", 2));

						G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*mW*src.get_val("VMIX", {ie+1, 0+1})*Gamma_UL[ae][ce]*src.get_val("WPARAM_SI_BSM", 10)-src.get_val("WPARAM_MATCH_BSM", {2,ce})*src.get_val("VMIX", {ie+1, 1+1})*Gamma_UR[ae][ce]*src.get_val("WPARAM_SI_BSM", 11))*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])/std::real(VCKM[2][1]));
					}

					if (ae < 3) {
						X_NL[ie][ae][be] = -g2 * src.get_val("VMIX", {ie+1, 0+1}) * Gamma_NL[ae][be]; 

						X_NR[ie][ae][be] = g2 * src.get_val("UMIX", {ie+1, 1+1}) * Gamma_NL[ae][be] * src.get_val("WPARAM_SI_BSM", {12, be}) / (sqrt(2.0) * mW * src.get_val("WPARAM_SI_BSM", 2)); //12 -> ME
					}

				}
			}
		}

		
		auto computeContributions = [&](int ie, auto func, double additionalFactor = 1.0) {
			double result = 0.0;
			for (int ae = 0; ae < 6; ++ae) {
				double msqOverMchSquared = std::pow(src.get_val("WPARAM_SI_BSM", {14, ae}) / src.get_val("WPARAM_SI_BSM", {13, ie}), 2.0);
				result += (X_UL[ie][ae][0] * X_UL[ie][ae][1] * func(msqOverMchSquared) + 
				src.get_val("WPARAM_SI_BSM", {13, ie}) / src.get_val("WPARAM_MATCH_SM", {5,1}) * X_UL[ie][ae][0] * X_UR[ie][ae][1] * func(msqOverMchSquared)) * additionalFactor;
			}
			return result;
		};


		auto hFunc10 = [](double x) { return h10(x); };
		auto hFunc20 = [](double x) { return h20(x); };
		auto hFunc50 = [](double x) { return h50(x); };
		auto hFunc60 = [](double x) { return h60(x); };

		double kappaFactor = -0.5 * src.get_val("WPARAM_SI_BSM", 6);

		
		for (int ie = 0; ie < 2; ++ie) {
			for (int je = 0; je < 2; ++je) {
				for (int ae = 0; ae < 6; ++ae) {
					double mchRatioSquared = pow(src.get_val("WPARAM_SI_BSM", {13, je}) / src.get_val("WPARAM_SI_BSM", {13, ie}), 2.0); //Mch : WPARAM_SI_BSM 13
					double msqOverMchSquared = pow(src.get_val("WPARAM_SI_BSM", {14, ae}) / src.get_val("WPARAM_SI_BSM", {13, ie}), 2.0);

					for (int be = 0; be < 3; ++be) {
						double msnOverMchSquared = pow(src.get_val("WPARAM_SI_BSM", {16, be}) / src.get_val("WPARAM_SI_BSM", {13, ie}), 2.0);
						B0c1 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (src.get_val("WPARAM_SI_BSM", {13, ie}) * src.get_val("WPARAM_SI_BSM", {13, ie})) * (0.5 * X_NL[ie][be][1] * X_NL[je][be][1] * f50(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));
						B0c2 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (src.get_val("WPARAM_SI_BSM", {13, ie}) * src.get_val("WPARAM_SI_BSM", {13, ie})) * (X_NR[ie][be][1] * X_NR[je][be][1] * std::fabs(src.get_val("WPARAM_SI_BSM", {13, je}) / src.get_val("WPARAM_SI_BSM", {13, ie})) * f60(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));	
					}

					C90c += X_UL[je][ae][1] * X_UL[ie][ae][2] * (2.0 * std::fabs(src.get_val("WPARAM_SI_BSM", {13, je}) / src.get_val("WPARAM_SI_BSM", {13, ie})) * f30(mchRatioSquared, msqOverMchSquared) * src.get_val("UMIX", {je+1, 0+1}) * src.get_val("UMIX", {ie+1,0+1}) - f40(mchRatioSquared, msqOverMchSquared) * src.get_val("VMIX", {je+1,0+1}) * src.get_val("VMIX", {ie+1,0+1}));

					if (ie == je)	{
						D90c += pow(mW / src.get_val("WPARAM_SI_BSM", {13, ie}), 2.0) * X_UL[ie][ae][1] * X_UL[ie][ae][2] * h30(msqOverMchSquared);
					}
				}
			}
		}

		for (int ie = 0; ie < 2; ++ie) {
			for (int ae = 0; ae < 6; ++ae) {
				for (int be = 0; be < 6; ++be) {
					double msqOverMchSquaredAe = pow(src.get_val("WPARAM_SI_BSM", {14, ae}) / src.get_val("WPARAM_SI_BSM", {13, ie}), 2.0);
					double msqOverMchSquaredBe = pow(src.get_val("WPARAM_SI_BSM", {14, be}) / src.get_val("WPARAM_SI_BSM", {13, ie}), 2.0);
					for (int ce = 0; ce < 3; ++ce) {
						C90c += X_UL[ie][be][1] * X_UL[ie][ae][2] * f40(msqOverMchSquaredAe, msqOverMchSquaredBe) * Gamma_UL[be][ce] * Gamma_UL[ae][ce];
					}
				}
			}
		}
		B90c = -(B0c1 - B0c2) * src.get_val("WPARAM_SI_BSM", 19) * std::pow(mW, 2.0) / (2.0 * std::pow(g2, 2.0));
		B100c = (B0c1 + B0c2) * src.get_val("WPARAM_SI_BSM", 19) * std::pow(mW, 2.0) / (2.0 * std::pow(g2, 2.0));
		C90c *= -src.get_val("WPARAM_SI_BSM", 19) / 8.0;
		D90c *= src.get_val("WPARAM_SI_BSM", 19);

		test = true;
		for (int ae = 0; ae < 6; ++ae) {
			if (!(std::fabs(src.get_val("WPARAM_SI_BSM", {14, ae})) > mW / 2. && std::fabs(src.get_val("WPARAM_SI_BSM", {15, ae})) > mW / 2.)) {
				test = false;
				break;
			}
		}


		for (int i = 0; i < Gamma_UL.size(); ++i) {
			for (int j = 0; j < Gamma_UL[i].size(); ++j) {
				dep_block->store_or_assign({1,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {1,i,j}}, Gamma_UL[i][j], 0., 0.)); //1
				dep_block->store_or_assign({2,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {2,i,j}}, Gamma_UR[i][j], 0., 0.)); //2
			}
		}

		for (int i = 0; i < X_UL.size(); ++i) {
			for (int j = 0; j < X_UL[i].size(); ++j) {
				for (int k = 0; k < X_UL[i][j].size(); ++k) {
					dep_block->store_or_assign({3,i,j,k}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {3,i,j,k}}, X_UL[i][j][k], 0., 0.)); //3
					dep_block->store_or_assign({4,i,j,k}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {4,i,j,k}}, X_UR[i][j][k], 0., 0.)); //4
					dep_block->store_or_assign({5,i,j,k}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {5,i,j,k}}, X_NL[i][j][k], 0., 0.)); //5
					dep_block->store_or_assign({6,i,j,k}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {6,i,j,k}}, X_NR[i][j][k], 0., 0.)); //6
				}
			}
		}

		for (int i = 0; i < Gamma_U.size(); ++i) {
			for (int j = 0; j < Gamma_U[i].size(); ++j) {
				dep_block->store_or_assign({7,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {7,i,j}}, Gamma_U[i][j], 0., 0.)); //7
				dep_block->store_or_assign({8,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {8,i,j}}, I_LR[i][j], 0., 0.)); //8
				dep_block->store_or_assign({9,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {9,i,j}}, P_U[i][j], 0., 0.)); //9
			}
		}

		for (int i = 0; i < Gamma_NL.size(); ++i) {
			for (int j = 0; j < Gamma_NL[i].size(); ++j) {
				dep_block->store_or_assign({10,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {10,i,j}}, Gamma_NL[i][j], 0., 0.)); //10
				dep_block->store_or_assign({11,i,j}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {11,i,j}}, Gamma_NR[i][j], 0., 0.)); //11
			}
		}

		for (int i = 0; i < G_aimn.size(); ++i) {
			for (int j = 0; j < G_aimn[i].size(); ++j) {
				for (int k = 0; k < G_aimn[i][j].size(); ++k) {
					for (int l = 0; l < G_aimn[i][j][k].size(); ++l) {
						dep_block->store_or_assign({12, i, j, k, l}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", {12, i, j, k, l}}, G_aimn[i][j][k][l], 0., 0.)); //12
					}
				}
			}
		}

        dep_block->store_or_assign(13, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", 13}, B90c, 0., 0.)); //13
		dep_block->store_or_assign(14, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", 14}, C90c, 0., 0.)); //14
		dep_block->store_or_assign(15, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", 15}, D90c, 0., 0.)); //15
		dep_block->store_or_assign(16, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "MATRIX_BSM", 16}, B100c, 0., 0.)); //16
    };

    iblock_c->compose_block("MATRIX_BSM", src_matrix, func_matrix);
}