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

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {

        src.at("ALPHA")->retrieve({})->get_val();

        double g2 = src.at("GAUGE")->retrieve(2)->get_val();
        double alpha_em = src.at("SMINPUTS")->retrieve(1)->get_val();

        double m_ds = src.at("MASS")->retrieve(1000001)->get_val();
        double m_us = src.at("MASS")->retrieve(1000002)->get_val();
        double m_ss = src.at("MASS")->retrieve(1000003)->get_val();
        double m_cs = src.at("MASS")->retrieve(1000004)->get_val();
        double m_bs = src.at("MASS")->retrieve(1000005)->get_val();
        double m_ts = src.at("MASS")->retrieve(1000006)->get_val();

        double m_d2s = src.at("MASS")->retrieve(2000001)->get_val();
        double m_u2s = src.at("MASS")->retrieve(2000002)->get_val();
        double m_s2s = src.at("MASS")->retrieve(2000003)->get_val();
        double m_c2s = src.at("MASS")->retrieve(2000004)->get_val();
        double m_b2s = src.at("MASS")->retrieve(2000005)->get_val();
        double m_t2s = src.at("MASS")->retrieve(2000006)->get_val();

        double m_gluino = src.at("MASS")->retrieve(1000021)->get_val();
        double m_c1 = src.at("MASS")->retrieve(1000024)->get_val();
        double m_c2 = src.at("MASS")->retrieve(1000037)->get_val();
        double mu_Q = src.at("HMIX")->retrieve(1)->get_val();

        double mqL3 = src.at("MSOFT")->retrieve(43)->get_val();
        double mbR = src.at("MSOFT")->retrieve(49)->get_val();

        double m_G = src.at("MASS")->retrieve(1000039)->get_val();

        std::map<int,int> neutralino = {{0, 1000022},{1, 1000023},{2, 1000025},{3, 1000035}};


        std::vector<double> m_neutralino = {src.at("MASS")->retrieve(neutralino[0])->get_val(), src.at("MASS")->retrieve(neutralino[1])->get_val(), src.at("MASS")->retrieve(neutralino[2])->get_val(), src.at("MASS")->retrieve(neutralino[3])->get_val()};
        double ad_22 = src.at("AD")->retrieve({3, 3})->get_val();
        double au_22 = src.at("AU")->retrieve({3, 3})->get_val();

        double yu_22 = src.at("YU")->retrieve({3, 3})->get_val();
        double yd_22 = src.at("YD")->retrieve({3, 3})->get_val();

        double sbot_mix_00 = src.at("SBOTMIX")->retrieve({0+1, 0+1})->get_val();
        double sbot_mix_01 = src.at("SBOTMIX")->retrieve({0+1, 1+1})->get_val();

        double stop_mix_00 = src.at("STOPMIX")->retrieve({0+1, 0+1})->get_val();
        double stop_mix_01 = src.at("STOPMIX")->retrieve({0+1, 1+1})->get_val();

        double umix_01 = src.at("UMIX")->retrieve({0+1, 1+1})->get_val();
        double umix_11 = src.at("UMIX")->retrieve({1+1, 1+1})->get_val();

        double vmix_01 = src.at("VMIX")->retrieve({0+1, 1+1})->get_val();
        double vmix_11 = src.at("VMIX")->retrieve({1+1, 1+1})->get_val();
        
        //0
        double sw2 = src.at("WPARAM_SI_SM")->retrieve(4)->get_val();

        // double alphas_MSOFT = QCDHelper::alpha_s(src.at("HMIX")->retrieve(0)->get_val()); // SUSY Breaking scale
        double alphas_MSOFT = QCDHelper::alpha_s(2.448e3); //TODO better
        double MSOFT = ParameterProxy(ParameterType::BSM).get_scale("MSOFT"); //TODO : better things to do with scale

        std::cout << "MSOFT" << MSOFT << std::endl; 
        double tan_beta = src.at("HMIX")->retrieve(2)->get_val();

        double factor = 2.0 / 3.0 * alphas_MSOFT / M_PI;

        double M_2 = src.at("MSOFT")->retrieve(2)->get_val();


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
                        src.at("NMIX")->retrieve({ie + 1, 3 + 1})->get_val() * src.at("NMIX")->retrieve({ie + 1, 2 + 1})->get_val() * 
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

        // return epsilonbp;
        
        //0p

        double epsilon0p = -2.0 / 3.0 * alphas_MSOFT / M_PI * 
                        (mu_Q + au_22 / tan_beta) / m_gluino *
                        (stop_mix_00 * stop_mix_00 * 
                            H2(m_t2s * m_t2s / m_gluino / m_gluino, m_ss * m_ss / m_gluino / m_gluino) +
                            stop_mix_01 * stop_mix_01 * 
                            H2(m_ts * m_ts / m_gluino / m_gluino, m_ss * m_ss / m_gluino / m_gluino));

        for(int ie = 0; ie < nb_neut; ++ie) {
            epsilon0p += yd_22 * yd_22 / 16.0 / M_PI / M_PI * 
                        src.at("NMIX")->retrieve({ie+1, 3+1})->get_val() * src.at("NMIX")->retrieve({ie+1, 2+1})->get_val() * 
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

        // return epsilon0p;

        //1p
        // Calcul du premier terme en utilisant yub[3], A_b, MqL3_Q, MbR_Q, mu_Q
        term1 = 1.0 / 16.0 / M_PI / M_PI * 
        (yd_22 * yd_22 * ad_22 / mu_Q * 
        H2(std::pow(mqL3 / mu_Q, 2), std::pow(mbR / mu_Q, 2))); //MbR_Q

        // Calcul du deuxième terme en utilisant g2, M2_Q, MqL3_Q, mu_Q
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

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double mW = src.at("MASS")->retrieve(24)->get_val();
		double alphas_mg = QCDHelper::alpha_s(src.at("MASS")->retrieve(1000021)->get_val());
		double ag = 1.0 - 7.0 / (12.0 * Pi) * alphas_mg;
		double aY = 1.0 + alphas_mg / (4.0 * Pi);

		// printf("ag : %.14lf\n",ag);
		// printf("aY : %.14lf\n",aY);

		// std::cout << "ag :" << ag << std::endl;
		// std::cout << "aY :" << aY << std::endl;
		// TODO : Ask Nazila Answer : keep complex
		double kappa = 1.0 / (pow(src.at("GAUGE")->retrieve(2)->get_val(), 2.) * 
						std::real((src.at("VCKM")->retrieve({2,2})->get_val())*(src.at("VCKM")->retrieve({2,1})->get_val()))); //VCKM 33 et 32

		// std::cout << "kappa :" << kappa << std::endl;
		// std::cout << "g2 :" << src.at("GAUGE")->retrieve(2)->get_val() << std::endl;
		// std::cout << "V22 :" << src.at("VCKM")->retrieve({2,2})->get_val() << std::endl;
		// std::cout << "V21 :" << src.at("VCKM")->retrieve({2,1})->get_val() << std::endl;	
		double kappaFactor = -0.5 * kappa;
		// std::cout << "g2 : " << src.at("GAUGE")->retrieve(2)->get_val() << std::endl;
		// std::cout << "kappa :" << kappa << std::endl;
		// std::cout << "aY :" << aY << std::endl;

		

		

		// std::cout << "VCKM 22 " << src.at("VCKM")->retrieve({2,2})->get_val() << std::endl;
		// std::cout << "VCKM 23 " << src.at("VCKM")->retrieve({2,1})->get_val() << std::endl;

		double tanb = src.at("HMIX")->retrieve(2)->get_val();
		double beta = std::atan(tanb);
		double z = pow(src.at("MASS")->retrieve(37)->get_val() / mW, 2.);
		double sinb = std::sin(std::atan(tanb));
		double cosb = std::cos(std::atan(tanb));
		double ct = src.at("STOPMIX")->retrieve({2,2})->get_val(); //TODO : 0 or 1 convention
		double st = src.at("STOPMIX")->retrieve({1,2})->get_val();  //TODO : 0 or 1 convention

		// std::cout << "cosb :" << cosb << std::endl;
		// std::cout << "mW :" << src.at("MASS")->retrieve(24)->get_val() << std::endl;

		// double ct = src.at("STOPMIX")->retrieve({1,1})->get_val(); //TODO : 0 or 1 convention
		// double st = src.at("STOPMIX")->retrieve({0,1})->get_val();  //TODO : 0 or 1 convention

        double lu = 1./tanb;
        double ld = -tanb;

		//TODO : IN progress
		Array1D_4 ME = {src.at("MASS")->retrieve(11)->get_val(), src.at("MASS")->retrieve(13)->get_val(), src.at("MASS")->retrieve(15)->get_val()}; 

		Array1D_3 Mch = {src.at("MASS")->retrieve(1000024)->get_val(), src.at("MASS")->retrieve(1000037)->get_val()};

		Array1D_7 MsqU = {src.at("MASS")->retrieve(1000002)->get_val(), src.at("MASS")->retrieve(1000004)->get_val(), src.at("MASS")->retrieve(1000006)->get_val(), 
				src.at("MASS")->retrieve(2000002)->get_val(), src.at("MASS")->retrieve(2000004)->get_val(), src.at("MASS")->retrieve(2000006)->get_val()};

		Array1D_7 MsqD = {src.at("MASS")->retrieve(1000001)->get_val(), src.at("MASS")->retrieve(1000003)->get_val(), src.at("MASS")->retrieve(1000005)->get_val(), 
				src.at("MASS")->retrieve(2000001)->get_val(), src.at("MASS")->retrieve(2000003)->get_val(), src.at("MASS")->retrieve(2000005)->get_val()};

		Array1D_4 Msn = {src.at("MASS")->retrieve(1000012)->get_val(), src.at("MASS")->retrieve(1000014)->get_val(), src.at("MASS")->retrieve(1000016)->get_val()};
		
		// std::cout << "mass 100002 : " <<  src.at("MASS")->retrieve(1000002)->get_val() << std::endl;
		// std::cout << "mass 100004 : " <<  src.at("MASS")->retrieve(1000004)->get_val() << std::endl;
		// for (int i  = 0; i<6; i++) {
		// 	std::cout << "MsqU[" << i+1 << "] = " << MsqU[i] << std::endl;
		// }
		// for(int ke=0;ke<6;ke++) printf("MsqU[%d] = %lf\n", ke,  MsqU[ke]);
		// for(int ke=0;ke<2;ke++) printf("Mch[%d] = %lf\n", ke,  Mch[ke]);
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
		// for(int ie=0;ie<3;ie++) printf("Msn[%d] : %.14lf\n",ie, Msn[ie]);
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

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {
        double yt= pow(src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val()/src.at("MASS")->retrieve(37)->get_val(),2.); // param->mass_H (25)
		Array1D_4 MU = {src.at("MASS")->retrieve(2)->get_val(), src.at("MASS")->retrieve(4)->get_val(), src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val()}; //TODO : size 3 not 4
		Array1D_4 MD = {src.at("MASS")->retrieve(2)->get_val(), src.at("MASS")->retrieve(3)->get_val(), src.at("WPARAM_MATCH_SM")->retrieve({5,1})->get_val()}; //TODO : size 3 not 4 // TODO : MD[0] -> mu like superiso but why ?

		// std::cout << "MU : " << src.at("MASS")->retrieve(2)->get_val() << std::endl;
		// std::cout << "MC : " << src.at("MASS")->retrieve(4)->get_val() << std::endl;
		// std::cout << "MT : " << src.at("WPARAM_MATCH_SM")->retrieve(6)->get_val() << std::endl;
		// std::cout << "MD : " << src.at("MASS")->retrieve(1)->get_val() << std::endl;
		// std::cout << "MS : " << src.at("MASS")->retrieve(3)->get_val() << std::endl;
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

    auto func_matrix = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {

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

		complex_t c11 = src.at("VCKM")->retrieve({0,0})->get_val();
		complex_t c12 = src.at("VCKM")->retrieve({0,1})->get_val();
		complex_t c13 = src.at("VCKM")->retrieve({0,2})->get_val();
		complex_t c21 = src.at("VCKM")->retrieve({1,0})->get_val();
		complex_t c22 = src.at("VCKM")->retrieve({1,1})->get_val();
		complex_t c23 = src.at("VCKM")->retrieve({1,2})->get_val();
		complex_t c31 = src.at("VCKM")->retrieve({2,0})->get_val();
		complex_t c32 = src.at("VCKM")->retrieve({2,1})->get_val();
		complex_t c33 = src.at("VCKM")->retrieve({2,2})->get_val();

		// printf("c11: %.14lf\n", c11.real());
		// printf("c12: %.14lf\n", c12.real());
		// printf("c13: %.14lf\n", c13.real());
		// printf("c21: %.14lf\n", c21.real());
		// printf("c22: %.14lf\n", c22.real());
		// printf("c23: %.14lf\n", c23.real());
		// printf("c31: %.14lf\n", c31.real());
		// printf("c32: %.14lf\n", c32.real());
		// printf("c33: %.14lf\n", c33.real());

		// complex_t c11 = src.at("RECKM")->retrieve(00)->get_val() + src.at("IMCKM")->retrieve(00)->get_val() * complex_t(0, 1);
		// complex_t c12 = src.at("RECKM")->retrieve(01)->get_val() + src.at("IMCKM")->retrieve(01)->get_val() * complex_t(0, 1);
		// complex_t c13 = src.at("RECKM")->retrieve(02)->get_val() + src.at("IMCKM")->retrieve(02)->get_val() * complex_t(0, 1);
		// complex_t c21 = src.at("RECKM")->retrieve(10)->get_val() + src.at("IMCKM")->retrieve(10)->get_val() * complex_t(0, 1);
		// complex_t c22 = src.at("RECKM")->retrieve(11)->get_val() + src.at("IMCKM")->retrieve(11)->get_val() * complex_t(0, 1);
		// complex_t c23 = src.at("RECKM")->retrieve(12)->get_val() + src.at("IMCKM")->retrieve(12)->get_val() * complex_t(0, 1);
		// complex_t c31 = src.at("RECKM")->retrieve(20)->get_val() + src.at("IMCKM")->retrieve(20)->get_val() * complex_t(0, 1);
		// complex_t c32 = src.at("RECKM")->retrieve(21)->get_val() + src.at("IMCKM")->retrieve(21)->get_val() * complex_t(0, 1);
		// complex_t c33 = src.at("RECKM")->retrieve(22)->get_val() + src.at("IMCKM")->retrieve(22)->get_val() * complex_t(0, 1);

		
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

        double mW = src.at("MASS")->retrieve(24)->get_val();
		double g2 = src.at("GAUGE")->retrieve(2)->get_val();
		// std::cout << "g2 : " << g2 << std::endl;
		if (src.at("WPARAM_SI_BSM")->retrieve(17)->get_val()) {
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
			Gamma_UL[2][2] = src.at("WPARAM_SI_BSM")->retrieve(4)->get_val();
			Gamma_UL[5][2] = -src.at("WPARAM_SI_BSM")->retrieve(5)->get_val();

			Gamma_UR[3][0] = 1.0;
			Gamma_UR[4][1] = 1.0;
			Gamma_UR[2][2] = src.at("WPARAM_SI_BSM")->retrieve(5)->get_val();
			Gamma_UR[5][2] = src.at("WPARAM_SI_BSM")->retrieve(4)->get_val();
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
				// printf("P_U[%d][%d] = %.8lf\n", ae, be, P_U[ae][be]);
			}
			
		}
		
		// std::cout << "UMIX[1][2] :" << src.at("UMIX")->retrieve({0+1, 1+1})->get_val() << std::endl;
		// std::cout << "UMIX[2][2] :" << src.at("UMIX")->retrieve({1+1, 1+1})->get_val() << std::endl;
		for (int ie = 0; ie < 2; ++ie) {
			for (int ae = 0; ae < 6; ++ae) {
				for (int be = 0; be < 3; ++be) {
					X_UL[ie][ae][be] = 0.0;
					X_UR[ie][ae][be] = 0.0;

					for (int ce = 0; ce < 3; ++ce) {
						X_UL[ie][ae][be] += -g2 * (
							src.at("WPARAM_SI_BSM")->retrieve(10)->get_val() * src.at("VMIX")->retrieve({ie+1, 0+1})->get_val() * Gamma_UL[ae][ce] -
							src.at("WPARAM_SI_BSM")->retrieve(11)->get_val() * src.at("VMIX")->retrieve({ie+1, 1+1})->get_val() * Gamma_UR[ae][ce] * src.at("WPARAM_MATCH_BSM")->retrieve({2,ce})->get_val() / (sqrt(2.0) * mW * src.at("WPARAM_SI_BSM")->retrieve(3)->get_val())
						) * std::real(VCKM[ce][be]);
						X_UR[ie][ae][be] += g2 * src.at("WPARAM_SI_BSM")->retrieve(11)->get_val() * src.at("UMIX")->retrieve({ie+1, 1+1})->get_val() * Gamma_UL[ae][ce] * std::real(VCKM[ce][be]) * src.at("WPARAM_MATCH_BSM")->retrieve({3,be})->get_val() / (sqrt(2.0) * mW * src.at("WPARAM_SI_BSM")->retrieve(2)->get_val());

						G_aimn[ae][ie][be][ce]=0.5/sqrt(2.)*(sqrt(2.)*mW*src.at("VMIX")->retrieve({ie+1, 0+1})->get_val()*Gamma_UL[ae][ce]*src.at("WPARAM_SI_BSM")->retrieve(10)->get_val()-src.at("WPARAM_MATCH_BSM")->retrieve({2,ce})->get_val()*src.at("VMIX")->retrieve({ie+1, 1+1})->get_val()*Gamma_UR[ae][ce]*src.at("WPARAM_SI_BSM")->retrieve(11)->get_val())*(std::real(VCKM[be][2])*std::real(VCKM[ce][1])/std::real(VCKM[2][2])/std::real(VCKM[2][1]));
					}

					if (ae < 3) {
						X_NL[ie][ae][be] = -g2 * src.at("VMIX")->retrieve({ie+1, 0+1})->get_val() * Gamma_NL[ae][be]; 

						X_NR[ie][ae][be] = g2 * src.at("UMIX")->retrieve({ie+1, 1+1})->get_val() * Gamma_NL[ae][be] * src.at("WPARAM_SI_BSM")->retrieve({12, be})->get_val() / (sqrt(2.0) * mW * src.at("WPARAM_SI_BSM")->retrieve(2)->get_val()); //12 -> ME
					}
					// printf("X_UR[%d][%d][%d] : %.14lf\n", ie+1, ae+1, be+1, X_UR[ie][ae][be]);
					// std::cout << "X_NR[" << ie+1 << "][" << ae+1 << "][" << be+1 << "] = " << X_NR[ie][ae][be] << std::endl;
					// std::cout << "X_UL[" << ie+1 << "][" << ae+1 << "][" << be+1 << "] = " << X_UL[ie][ae][be] << std::endl;
				}
			}
		}
		// printf("param->g2 : %.14lf\n", g2);
		// printf("param->g2 : %.14lf\n", g2);
		// printf("param->mW : %.14lf\n", mW);
		// printf("cosb : %.14lf\n", src.at("WPARAM_SI_BSM")->retrieve(2)->get_val().real());
		// for (int i = 0; i < 6; ++i) {
		// 	for (int j = 0; j < 3; ++j) {
		// 		printf("Gamma_UL[%d][%d] = %f\n", i+1, j+1, Gamma_UL[i][j]);
		// 	}
		// }

		// for (int i = 0; i<3; i++) {
		// 	for (int j = 0; j<3; j++) {
		// 		printf("VCKM[%d][%d] = %f\n",  i, j, std::real(VCKM[i][j]));
		// 	}
		// }
		
		auto computeContributions = [&](int ie, auto func, double additionalFactor = 1.0) {
			double result = 0.0;
			for (int ae = 0; ae < 6; ++ae) {
				double msqOverMchSquared = std::pow(src.at("WPARAM_SI_BSM")->retrieve({14, ae})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);
				result += (X_UL[ie][ae][0] * X_UL[ie][ae][1] * func(msqOverMchSquared) + 
				src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val() / src.at("WPARAM_MATCH_SM")->retrieve({5,1})->get_val() * X_UL[ie][ae][0] * X_UR[ie][ae][1] * func(msqOverMchSquared)) * additionalFactor;
			}
			return result;
		};


		auto hFunc10 = [](double x) { return h10(x); };
		auto hFunc20 = [](double x) { return h20(x); };
		auto hFunc50 = [](double x) { return h50(x); };
		auto hFunc60 = [](double x) { return h60(x); };

		double kappaFactor = -0.5 * src.at("WPARAM_SI_BSM")->retrieve(6)->get_val();

		
		for (int ie = 0; ie < 2; ++ie) {
			for (int je = 0; je < 2; ++je) {
				for (int ae = 0; ae < 6; ++ae) {
					double mchRatioSquared = pow(src.at("WPARAM_SI_BSM")->retrieve({13, je})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0); //Mch : WPARAM_SI_BSM 13
					double msqOverMchSquared = pow(src.at("WPARAM_SI_BSM")->retrieve({14, ae})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);

					for (int be = 0; be < 3; ++be) {
						double msnOverMchSquared = pow(src.at("WPARAM_SI_BSM")->retrieve({16, be})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);
						B0c1 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val() * src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val()) * (0.5 * X_NL[ie][be][1] * X_NL[je][be][1] * f50(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));
						B0c2 += X_UL[je][ae][1] * X_UL[ie][ae][2] / (src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val() * src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val()) * (X_NR[ie][be][1] * X_NR[je][be][1] * std::fabs(src.at("WPARAM_SI_BSM")->retrieve({13, je})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val()) * f60(mchRatioSquared, msqOverMchSquared, msnOverMchSquared));	
					}

					C90c += X_UL[je][ae][1] * X_UL[ie][ae][2] * (2.0 * std::fabs(src.at("WPARAM_SI_BSM")->retrieve({13, je})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val()) * f30(mchRatioSquared, msqOverMchSquared) * src.at("UMIX")->retrieve({je+1, 0+1})->get_val() * src.at("UMIX")->retrieve({ie+1,0+1})->get_val() - f40(mchRatioSquared, msqOverMchSquared) * src.at("VMIX")->retrieve({je+1,0+1})->get_val() * src.at("VMIX")->retrieve({ie+1,0+1})->get_val());

					if (ie == je)	{
						D90c += pow(mW / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0) * X_UL[ie][ae][1] * X_UL[ie][ae][2] * h30(msqOverMchSquared);
					}
				}
			}
		}

		for (int ie = 0; ie < 2; ++ie) {
			for (int ae = 0; ae < 6; ++ae) {
				for (int be = 0; be < 6; ++be) {
					double msqOverMchSquaredAe = pow(src.at("WPARAM_SI_BSM")->retrieve({14, ae})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);
					double msqOverMchSquaredBe = pow(src.at("WPARAM_SI_BSM")->retrieve({14, be})->get_val() / src.at("WPARAM_SI_BSM")->retrieve({13, ie})->get_val(), 2.0);
					for (int ce = 0; ce < 3; ++ce) {
						C90c += X_UL[ie][be][1] * X_UL[ie][ae][2] * f40(msqOverMchSquaredAe, msqOverMchSquaredBe) * Gamma_UL[be][ce] * Gamma_UL[ae][ce];
					}
				}
			}
		}
		B90c = -(B0c1 - B0c2) * src.at("WPARAM_SI_BSM")->retrieve(19)->get_val() * std::pow(mW, 2.0) / (2.0 * std::pow(g2, 2.0));
		B100c = (B0c1 + B0c2) * src.at("WPARAM_SI_BSM")->retrieve(19)->get_val() * std::pow(mW, 2.0) / (2.0 * std::pow(g2, 2.0));
		C90c *= -src.at("WPARAM_SI_BSM")->retrieve(19)->get_val() / 8.0;
		D90c *= src.at("WPARAM_SI_BSM")->retrieve(19)->get_val();

		test = true;
		for (int ae = 0; ae < 6; ++ae) {
			if (!(std::fabs(src.at("WPARAM_SI_BSM")->retrieve({14, ae})->get_val()) > mW / 2. && std::fabs(src.at("WPARAM_SI_BSM")->retrieve({15, ae})->get_val()) > mW / 2.)) {
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