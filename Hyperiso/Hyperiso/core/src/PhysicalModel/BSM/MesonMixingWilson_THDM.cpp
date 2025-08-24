#include "MesonMixingWilson_THDM.h"

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();;
	// if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	// else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta=src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	
	double D0h,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
	}

}

C_mix_bd_1_tilde_THDM::C_mix_bd_1_tilde_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_tilde_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	// if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	// else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	
	scalar_t Cp1_chargedhiggs=0.;

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));

	}

}

C_mix_bd_2_THDM::C_mix_bd_2_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_2_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	// if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	// else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C2_chargedhiggs=0.;

	
	double D0h,D0h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}

C_mix_bd_1_THDM::C_mix_bd_1_THDM() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
            {ParameterType::SM, "MASS", 3},                       //m_s
            {ParameterType::SM, "GAUGE", 2},                       //g_2 
            {ParameterType::SM, "VCKM", LhaID(0, 0)},             // V_ud
            {ParameterType::SM, "VCKM", LhaID(0, 1)},             // V_us
            {ParameterType::SM, "VCKM", LhaID(0, 2)},             // V_ub
            {ParameterType::SM, "VCKM", LhaID(1, 0)},             // V_cd
            {ParameterType::SM, "VCKM", LhaID(1, 1)},             // V_cs
            {ParameterType::SM, "VCKM", LhaID(1, 2)},             // V_cb
            {ParameterType::SM, "VCKM", LhaID(2, 0)},             // V_td
            {ParameterType::SM, "VCKM", LhaID(2, 1)},             // V_ts
            {ParameterType::SM, "VCKM", LhaID(2, 2)},             // V_tb
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 0)
    };
}

double C_mix_bd_1_THDM::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q;
	if(gen==1) m_q=src.at({ParameterType::SM, "MASS", 1})->get_val();
	else m_q=src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta= src.at({ParameterType::BSM, "MINPAR", 3})->get_val(); //tan_b
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
	// m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	// m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	scalar_t C2_chargedhiggs=0.;
	scalar_t C3_chargedhiggs=0.;
	scalar_t C4_chargedhiggs=0.;
	scalar_t C5_chargedhiggs=0.;
	
	scalar_t Cp1_chargedhiggs=0.;
	scalar_t Cp2_chargedhiggs=0.;
	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][gen])*conj(V_CKM[j][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}

}