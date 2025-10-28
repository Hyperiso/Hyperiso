#include "MesonMixingWilsonSUSY.h"

C_mix_bd_1_SUSY::C_mix_bd_1_SUSY() : WilsonCoefficient("C_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(1050105, 4141, 0, 1)
    };
}

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added 0eration dependence, V_CKM[ie][2] -> V_CKM[ie][0] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C1_gluino=0.;

	double D2g,D0g;

	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie] and Z_D[5][ie] -> Z_D[0+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[2][i]*Z_D[2][j]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[0][i])*conj(Z_D[0][j])/(144.*pow(PI, 2));
		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	
	scalar_t C1_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added 0eration dependence, Yd[2] -> Yd[0] and V_CKM[Ke][2] -> V_CKM[Ke][0] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   
						
						C1_chargino += -D2ch*pow(V_CKM[k][2], 2)*(-Q_e*Z_p[0][a]*conj(Z_U[k][i])*swi + Yu[k]*Z_p[1][a]*conj(Z_U[k+3][i]))*(-Q_e*Z_p[0][b]*conj(Z_U[k][j])*swi + Yu[k]*Z_p[1][b]*conj(Z_U[k+3][j]))*(-Q_e*Z_U[k][i]*conj(Z_p[0][b])*swi + Z_U[k+3][i]*conj(Yu[k])*conj(Z_p[1][b]))*(-Q_e*Z_U[k][j]*conj(Z_p[0][a])*swi + Z_U[k+3][j]*conj(Yu[k])*conj(Z_p[1][a]))*pow(conj(V_CKM[k][0]), 2)/(32.0*pow(PI, 2));
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C1_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					C1_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*
					(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][b])*
					(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*
					(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][b]))/(64.0*pow(PI, 2)) - 
					D2ne*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*
					(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][b])*
					(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*
					(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][b]))/(32.0*pow(PI, 2));
				}
			}
		}
	}
	


    //mixed ->


	scalar_t C1_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 
				
				C1_mixed += -D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[0][i]*Z_D[0][j]*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[2][i])/(2.0*cw*sw) + 
				conj(Yd[2])*conj(Z_D[5][i])*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[2][j])/(2.0*cw*sw) + conj(Yd[2])*conj(Z_D[5][j])*conj(Z_N[2][a])) + 
				Z_D[2][i]*Z_D[2][j]*pow(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]), 2))/(96.0*pow(PI, 2)) - 
				D2mix*Z_D[2][j]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + 
				conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*conj(Z_D[0][i])/(8.0*pow(PI, 2));
    		}
		}
	}

    //higgs PInguin ->




	return C1_chargedhiggs+C1_gluino+C1_chargino+C1_neutralino;

}

//BD_1_tilde

C_mix_bd_1_tilde_SUSY::C_mix_bd_1_tilde_SUSY() : WilsonCoefficient("CT_BD_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(1050105, 4242, 0, 1)
    };
}

double C_mix_bd_1_tilde_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	
	scalar_t Cp1_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added 0eration dependence, V_CKM[ie][2] -> V_CKM[ie][0] */
        
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t Cp1_gluino=0.;

	double D2g,D0g;

	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie] and Z_D[5][ie] -> Z_D[0+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		Cp1_gluino += -Z_D[5][i]*Z_D[5][j]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[0+3][i])*conj(Z_D[0+3][j])/(144.*pow(PI, 2.));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	
	scalar_t Cp1_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added 0eration dependence, Yd[2] -> Yd[0] and V_CKM[Ke][2] -> V_CKM[Ke][0] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
			
						Cp1_chargino += -D2ch*pow(V_CKM[k][2], 2)*pow(Yd[2], 2)*Z_m[1][a]*Z_m[1][b]*Z_U[k][i]*Z_U[k][j]*pow(conj(V_CKM[k][0]), 2)*pow(conj(Yd[0]), 2)*conj(Z_m[1][a])*conj(Z_m[1][b])*conj(Z_U[k][i])*conj(Z_U[k][j])/(32.0*pow(PI, 2));

					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t Cp1_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					Cp1_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][b])/(3.0*cw) + 
					Yd[2]*Z_D[2][j]*conj(Z_N[2][b]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[0+3][i])/(3.0*cw) + 
					Z_N[2][b]*conj(Yd[0])*conj(Z_D[0][i]))/(64.0*pow(PI, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][b])/(3.0*cw) + Yd[2]*Z_D[2][j]*conj(Z_N[2][b]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + 
					Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[0+3][i])/(3.0*cw) + Z_N[2][b]*conj(Yd[0])*conj(Z_D[0][i]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][i]*Z_N[2][a])/(32.0*pow(PI, 2));

				}
			}
		}
	}
	


    //mixed ->


	scalar_t Cp1_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 	

				Cp1_mixed += -D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[0+3][i]*Z_D[0+3][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[5][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[2])*conj(Z_D[2][i]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[5][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[2])*conj(Z_D[2][j])) + Z_D[5][i]*Z_D[5][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][i]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j])))/(96.0*pow(PI, 2)) - D2mix*Z_D[5][j]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*conj(Z_D[0+3][i])/(8.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->



}

//BD2

C_mix_bd_2_SUSY::C_mix_bd_2_SUSY() : WilsonCoefficient("C_BD_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(1050105, 3131, 0, 1)
    };
}

double C_mix_bd_2_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }

	scalar_t C2_chargedhiggs=0.;

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added 0eration dependence, V_CKM[ie][2] -> V_CKM[ie][0] */
        
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C2_gluino=0.;

	double D2g,D0g;

	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie] and Z_D[5][ie] -> Z_D[0+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[2][j]*pow(g_3, 4.)*conj(Z_D[0+3][i])*conj(Z_D[0+3][j])/(288.*pow(PI, 2));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	
	scalar_t C2_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 0eration dependence, Yd[2] -> Yd[0] and V_CKM[Ke][2] -> V_CKM[Ke][0] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   
						C2_chargino += 0.;

					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C2_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					C2_neutralino += D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[0+3][i])/(3.0*cw) + Z_N[2][b]*conj(Yd[0])*conj(Z_D[0][i]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[0+3][j])/(3.0*cw) + 
					Z_N[2][b]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/
					(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][a])/(32.0*pow(PI, 2));
					
				}
			}
		}
	}
	


    //mixed ->

	scalar_t C2_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 
				

				C2_mixed += D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[2][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][i]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*conj(Z_D[2][j]) + 3.0*Z_D[2][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - 
					Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[0+3][i]))/(48.0*pow(PI, 2)) + D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][a])*conj(Z_D[0+3][i])*conj(Z_D[0+3][j])/(48.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->

}

//BD2_tild

C_mix_bd_2_tilde_SUSY::C_mix_bd_2_tilde_SUSY() : WilsonCoefficient("CT_BD_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(1050105, 3232, 0, 1)
    };
}

double C_mix_bd_2_tilde_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  

	scalar_t Cp2_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added 0eration dependence, V_CKM[ie][2] -> V_CKM[ie][0] */
        
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));
	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t Cp2_gluino=0.;

	double D2g,D0g;

	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie] and Z_D[5][ie] -> Z_D[0+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[5][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[0][i])*conj(Z_D[0][j])/(288.*pow(PI, 2.));
		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	


	scalar_t Cp2_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 0eration dependence, Yd[2] -> Yd[0] and V_CKM[Ke][2] -> V_CKM[Ke][0] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   
						
						Cp2_chargino += 0.;
						
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C1_neutralino=0.;
	scalar_t C2_neutralino=0.;
	scalar_t C3_neutralino=0.;
	scalar_t C4_neutralino=0.;
	scalar_t C5_neutralino=0.;
	scalar_t Cp1_neutralino=0.;
	scalar_t Cp2_neutralino=0.;
	scalar_t Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					C1_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*
					(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][b])*
					(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*
					(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][b]))/(64.0*pow(PI, 2)) - 
					D2ne*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*
					(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][b])*
					(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*
					(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][b]))/(32.0*pow(PI, 2));
					

					Cp2_neutralino += D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
					Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][b]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - 
					conj(Z_N[1][b])*cw)*conj(Z_D[0][j])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][j])*conj(Z_N[2][b]))/(32.0*pow(PI, 2));

				}
			}
		}
	}
	


    //mixed ->


	scalar_t Cp2_mixed=0.;


	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 

				Cp2_mixed += D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[5][i]*pow(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + 
				conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]), 2)*conj(Z_D[5][j]) + 3.0*Z_D[5][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*conj(Z_D[0][i]))/(48.0*pow(PI, 2)) + 
				D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*conj(Z_D[0][i])*conj(Z_D[0][j])/(48.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->


}

//BD3

C_mix_bd_3_SUSY::C_mix_bd_3_SUSY() : WilsonCoefficient("C_BD_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(1050105, 7171, 0, 1)
    };
}

double C_mix_bd_3_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  

	scalar_t C3_chargedhiggs=0.;

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added 0eration dependence, V_CKM[ie][2] -> V_CKM[ie][0] */
        
		C3_chargedhiggs += 0.;

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C3_gluino=0.;

	double D2g,D0g;

	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie] and Z_D[5][ie] -> Z_D[0+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		C3_gluino += -D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[2][j]*pow(g_3,4.)*conj(Z_D[0+3][i])*conj(Z_D[0+3][j])/(96.*pow(PI, 2));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	scalar_t C3_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 0eration dependence, Yd[2] -> Yd[0] and V_CKM[Ke][2] -> V_CKM[Ke][0] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   
						C3_chargino += -D0ch*M_ch[a]*M_ch[b]*pow(V_CKM[k][2], 2)*Z_m[1][a]*Z_m[1][b]*Z_U[k][i]*Z_U[k][j]*(-Q_e*Z_p[0][a]*conj(Z_U[k][i])*swi + Yu[k]*Z_p[1][a]*conj(Z_U[k+3][i]))*(-Q_e*Z_p[0][b]*conj(Z_U[k][j])*swi + Yu[k]*Z_p[1][b]*conj(Z_U[k+3][j]))*pow(conj(V_CKM[k][0]), 2)*pow(conj(Yd[0]), 2)/(32.0*pow(PI, 2));

					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C3_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					
					C3_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*((-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][j]*Z_N[2][b]) - (-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][b]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][j]*Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[0+3][i])/(3.0*cw) + Z_N[2][b]*conj(Yd[0])*conj(Z_D[0][i]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][i]*Z_N[2][a])/(32.0*pow(PI, 2));

				}
			}
		}
	}
	


    //mixed ->


	scalar_t C3_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 
				
				C3_mixed += D0mix*M_ch0[a]*M_g*Z_D[2][i]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][i]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*conj(Z_D[2][j])/(48.0*pow(PI, 2)) - D0mix*M_ch0[a]*M_g*Z_D[2][j]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[0+3][i])/(48.0*pow(PI, 2)) + 
				D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
				Yd[2]*Z_D[5][j]*Z_N[2][a])*conj(Z_D[0+3][i])*conj(Z_D[0+3][j])/(48.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->

}

//BD3_tilde

C_mix_bd_3_tilde_SUSY::C_mix_bd_3_tilde_SUSY() : WilsonCoefficient("CT_BD_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(1050105, 7272, 0, 1)
    };
}

double C_mix_bd_3_tilde_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  

	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added 0eration dependence, V_CKM[ie][2] -> V_CKM[ie][0] */
        
		Cp3_chargedhiggs += 0.;
	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	scalar_t Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie] and Z_D[5][ie] -> Z_D[0+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[5][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[0][i])*conj(Z_D[0][j])/(96.*pow(PI, 2));
		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	scalar_t Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added 0eration dependence, Yd[2] -> Yd[0] and V_CKM[Ke][2] -> V_CKM[Ke][0] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						
						Cp3_chargino += -D0ch*M_ch[a]*M_ch[b]*pow(V_CKM[k][2], 2)*pow(Yd[2], 2)*(-Q_e*Z_U[k][i]*conj(Z_p[0][b])*swi + Z_U[k+3][i]*conj(Yu[k])*conj(Z_p[1][b]))*(-Q_e*Z_U[k][j]*conj(Z_p[0][a])*swi + Z_U[k+3][j]*conj(Yu[k])*conj(Z_p[1][a]))*pow(conj(V_CKM[k][0]), 2)*conj(Z_m[1][a])*conj(Z_m[1][b])*conj(Z_U[k][i])*conj(Z_U[k][j])/(32.0*pow(PI, 2));
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					Cp3_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*(-(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[0][j])/(2.0*cw*sw) + 
					conj(Yd[0])*conj(Z_D[0+3][j])*conj(Z_N[2][b])) + (-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][b])/(3.0*cw) + Yd[2]*Z_D[2][j]*conj(Z_N[2][b]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
					conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a])))*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + 
					Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][b]))/(32.0*pow(PI, 2));

				}
			}
		}
	}
	


    //mixed ->


	scalar_t Cp3_mixed=0.;


	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 
				
				Cp3_mixed += D0mix*M_ch0[a]*M_g*Z_D[5][i]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + 
				conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]), 2)*conj(Z_D[5][j])/(48.0*pow(PI, 2)) - D0mix*M_ch0[a]*M_g*Z_D[5][j]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + 
				conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*conj(Z_D[0][i])/(48.0*pow(PI, 2)) + D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*conj(Z_D[0][i])*conj(Z_D[0][j])/(48.0*pow(PI, 2));
    		}
		}
	}

    //higgs PInguin ->

}

//BD4

C_mix_bd_4_SUSY::C_mix_bd_4_SUSY() : WilsonCoefficient("C_BD_4", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(1050105, 3132, 0, 1)
    };
}

double C_mix_bd_4_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C4_chargedhiggs=0.;

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added 0eration dependence, V_CKM[ie][2] -> V_CKM[ie][0] */
        
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C4_gluino=0.;

	double D2g,D0g;

	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie] and Z_D[5][ie] -> Z_D[0+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[0][i])*conj(Z_D[0+3][j])/(48.*pow(PI, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[0][i])*conj(Z_D[0+3][j]) + 11.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[0][j])*conj(Z_D[0+3][i]))/(72.*pow(PI, 2));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	
	scalar_t C4_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 0eration dependence, Yd[2] -> Yd[0] and V_CKM[Ke][2] -> V_CKM[Ke][0] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   
						C4_chargino += D2ch*pow(V_CKM[k][2], 2)*Yd[2]*Z_m[1][a]*Z_U[k][j]*(-Q_e*Z_p[0][b]*conj(Z_U[k][j])*swi + Yu[k]*Z_p[1][b]*conj(Z_U[k+3][j]))*(-Q_e*Z_U[k][i]*conj(Z_p[0][b])*swi + Z_U[k+3][i]*conj(Yu[k])*conj(Z_p[1][b]))*pow(conj(V_CKM[k][0]), 2)*conj(Yd[0])*conj(Z_m[1][a])*conj(Z_U[k][i])/(8.0*pow(PI, 2));

					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C4_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][j]*Z_N[2][b]) + (-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][b]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][j]*Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + 
					conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][b]))/(8.0*pow(PI, 2));


				}
			}
		}
	}
	


    //mixed ->

	scalar_t C4_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 
				
				C4_mixed += D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[2][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*conj(Z_D[0+3][i]) + Z_D[5][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[0][i]))/(16.0*pow(PI, 2)) - 
				D2mix*pow(g_3, 2)*(-3.0*Z_D[2][j]*Z_D[5][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][i]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a])) - 3.0*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[5][i])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[2])*conj(Z_D[2][i]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][a])*conj(Z_D[0][j])*conj(Z_D[0+3][i]))/(24.0*pow(PI, 2)) - 
				D2mix*pow(g_3, 2)*(-Z_D[2][j]*Z_D[5][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a])) - (-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[0][j])*conj(Z_D[0+3][i]))/(24.0*pow(PI, 2)) - 
				D2mix*pow(g_3, 2)*(Z_D[2][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*conj(Z_D[0][i]) + Z_D[5][j]*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
				Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*conj(Z_D[0+3][i]))/(24.0*pow(PI, 2));


    		}
		}
	}

    //higgs PInguin ->




	scalar_t delta_d[6][6],delta_d_LL[3][3],delta_d_LR[3][3],delta_d_RL[3][3],delta_d_RR[3][3];
	scalar_t delta_u[6][6],delta_u_LL[3][3],delta_u_LR[3][3],delta_u_RL[3][3],delta_u_RR[3][3];
	
	double m_av = (M_U[0]+M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_D[0]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double M_A=src.at({ParameterType::BSM, "MASS", 36})->get_val();
	double A_t=src.at({ParameterType::BSM, "AU", LhaID(3,3)})->get_val(); //A_t
	double mu = src.at({ParameterType::BSM, "HMIX", 1})->get_val();
	scalar_t M_2 = src.at({ParameterType::BSM, "MSOFT", 2})->get_val();
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	scalar_t x_2 = pow(abs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*PI) ;
	double alpha_2 = pow(g_2,2.)/(4.*PI);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*PI*pow(m_av,2));
	
	
	scalar_t V_tb = V_CKM[2][2]; 
	scalar_t V_tq = V_CKM[2][0];
	scalar_t a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*PI*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	scalar_t a2 = -alpha_s*pow(M_g,2)*delta_d_LL[2][0]*delta_d_RR[2][0]*pow(h1(x_g),2);
	scalar_t a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[2][0]*V_tb*conj(V_tq)/(pow(M_W,2));
	scalar_t a4 = alpha_2*M_2*M_g*delta_u_LL[2][0]*delta_d_RR[2][0]*h1(x_g)*h4(x_2,x_g);

	scalar_t C4_higgspenguin = a1*(a2+a3+a4);

}

//BD_5

C_mix_bd_5_SUSY::C_mix_bd_5_SUSY() : WilsonCoefficient("C_BD_5", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 1},                       //m_d 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(1050105, 7172, 0, 1)
    };
}

double C_mix_bd_5_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 1})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C5_chargedhiggs=0.;
	

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][0])*conj(V_CKM[j][0]); 	/* NM: added 0eration dependence, V_CKM[ie][2] -> V_CKM[ie][0] */
        
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C5_gluino=0.;

	double D2g,D0g;

	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie] and Z_D[5][ie] -> Z_D[0+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		C5_gluino += -D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[0][i])*conj(Z_D[0+3][j])/(144.*pow(PI, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[0][i])*conj(Z_D[0+3][j]) + 3.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[0][j])*conj(Z_D[0+3][i]))/(72.*pow(PI, 2));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	


	scalar_t C5_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 0eration dependence, Yd[2] -> Yd[0] and V_CKM[Ke][2] -> V_CKM[Ke][0] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   
						C5_chargino += -D0ch*M_ch[a]*M_ch[b]*pow(V_CKM[k][2], 2)*Yd[2]*Z_m[1][b]*Z_U[k][i]*(-Q_e*Z_p[0][b]*conj(Z_U[k][j])*swi + Yu[k]*Z_p[1][b]*conj(Z_U[k+3][j]))*(-Q_e*Z_U[k][j]*conj(Z_p[0][a])*swi + Z_U[k+3][j]*conj(Yu[k])*conj(Z_p[1][a]))*pow(conj(V_CKM[k][0]), 2)*conj(Yd[0])*conj(Z_m[1][a])*conj(Z_U[k][i])/(16.0*pow(PI, 2));

					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C5_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					C5_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[0+3][i])/(3.0*cw) + 
					Z_N[2][b]*conj(Yd[0])*conj(Z_D[0][i]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][b])*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
					conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))/(16.0*pow(PI, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
					Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][b]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0- Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][b]))/(8.0*pow(PI, 2));

				}
			}
		}
	}
	


    //mixed ->


	scalar_t C5_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 0eration dependence, Z_D[2][ie] -> Z_D[0][ie], Z_D[5][ie] -> Z_D[0+3][ie] and Yd[2] -> Yd[0] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 
				
				C5_mixed += -D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[2][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*conj(Z_D[0+3][i]) + Z_D[5][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[0][i]))/(48.0*pow(PI, 2)) + 
				D2mix*pow(g_3, 2)*(-Z_D[2][j]*Z_D[5][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][i]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a])) - (-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[5][i])/(3*cw) + 
				Z_N[2][a]*conj(Yd[2])*conj(Z_D[2][i]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][a])*conj(Z_D[0][j])*conj(Z_D[0+3][i]))/(24.0*pow(PI, 2)) + 
				D2mix*pow(g_3, 2)*(-Z_D[2][j]*Z_D[5][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a])) - (-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[0][j])*conj(Z_D[0+3][i]))/(8.0*pow(PI, 2)) + 
				D2mix*pow(g_3, 2)*(Z_D[2][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[0+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[0])*conj(Z_D[0][j]))*conj(Z_D[0][i]) + Z_D[5][j]*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[0][i])/(2.0*cw*sw) + conj(Yd[0])*conj(Z_D[0+3][i])*conj(Z_N[2][a]))*conj(Z_D[0+3][i]))/(8.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->

}

//BS


C_mix_bs_1_SUSY::C_mix_bs_1_SUSY() : WilsonCoefficient("C_BS_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 3},                       //m_s 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(3050305, 4141, 0, 1)
    };
}

double C_mix_bs_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C1_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][1])*conj(V_CKM[j][1]); 	/* NM: added 1eration dependence, V_CKM[ie][2] -> V_CKM[ie][1] */
        
		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[i]*m_u_pow_2[j]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(PI,2.)*pow(M_W,4.));
	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C1_gluino=0.;

	double D2g,D0g;

	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie] and Z_D[5][ie] -> Z_D[1+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[2][i]*Z_D[2][j]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[1][i])*conj(Z_D[1][j])/(144.*pow(PI, 2));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	
	scalar_t C1_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 1eration dependence, Yd[2] -> Yd[1] and V_CKM[Ke][2] -> V_CKM[Ke][1] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   
						
						C1_chargino += -D2ch*pow(V_CKM[k][2], 2)*(-Q_e*Z_p[0][a]*conj(Z_U[k][i])*swi + Yu[k]*Z_p[1][a]*conj(Z_U[k+3][i]))*(-Q_e*Z_p[0][b]*conj(Z_U[k][j])*swi + Yu[k]*Z_p[1][b]*conj(Z_U[k+3][j]))*(-Q_e*Z_U[k][i]*conj(Z_p[0][b])*swi + Z_U[k+3][i]*conj(Yu[k])*conj(Z_p[1][b]))*(-Q_e*Z_U[k][j]*conj(Z_p[0][a])*swi + Z_U[k+3][j]*conj(Yu[k])*conj(Z_p[1][a]))*pow(conj(V_CKM[k][1]), 2)/(32.0*pow(PI, 2));
						
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C1_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					C1_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*
					(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][b])*
					(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))*
					(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][b]))/(64.0*pow(PI, 2)) - 
					D2ne*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*
					(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][b])*
					(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))*
					(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][b]))/(32.0*pow(PI, 2));

				}
			}
		}
	}
	


    //mixed ->


	scalar_t C1_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 
				
				C1_mixed += -D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[1][i]*Z_D[1][j]*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[2][i])/(2.0*cw*sw) + 
				conj(Yd[2])*conj(Z_D[5][i])*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[2][j])/(2.0*cw*sw) + conj(Yd[2])*conj(Z_D[5][j])*conj(Z_N[2][a])) + 
				Z_D[2][i]*Z_D[2][j]*pow(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]), 2))/(96.0*pow(PI, 2)) - 
				D2mix*Z_D[2][j]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + 
				conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))*conj(Z_D[1][i])/(8.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->

}

//BS_1_tilde 

C_mix_bs_1_tilde_SUSY::C_mix_bs_1_tilde_SUSY() : WilsonCoefficient("CT_BS_1", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 3},                       //m_s 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(3050305, 4242, 0, 1)
    };
}

double C_mix_bs_1_tilde_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	
	scalar_t Cp1_chargedhiggs=0.;

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][1])*conj(V_CKM[j][1]); 	/* NM: added 1eration dependence, V_CKM[ie][2] -> V_CKM[ie][1] */
        
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(PI,2.)*pow(M_W,4.));
		
	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t Cp1_gluino=0.;

	double D2g,D0g;

	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie] and Z_D[5][ie] -> Z_D[1+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		Cp1_gluino += -Z_D[5][i]*Z_D[5][j]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[1+3][i])*conj(Z_D[1+3][j])/(144.*pow(PI, 2.));
		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	
	scalar_t Cp1_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 1eration dependence, Yd[2] -> Yd[1] and V_CKM[Ke][2] -> V_CKM[Ke][1] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   
						Cp1_chargino += -D2ch*pow(V_CKM[k][2], 2)*pow(Yd[2], 2)*Z_m[1][a]*Z_m[1][b]*Z_U[k][i]*Z_U[k][j]*pow(conj(V_CKM[k][1]), 2)*pow(conj(Yd[1]), 2)*conj(Z_m[1][a])*conj(Z_m[1][b])*conj(Z_U[k][i])*conj(Z_U[k][j])/(32.0*pow(PI, 2));
						

						
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t Cp1_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					Cp1_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][b])/(3.0*cw) + 
					Yd[2]*Z_D[2][j]*conj(Z_N[2][b]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[1+3][i])/(3.0*cw) + 
					Z_N[2][b]*conj(Yd[1])*conj(Z_D[1][i]))/(64.0*pow(PI, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][b])/(3.0*cw) + Yd[2]*Z_D[2][j]*conj(Z_N[2][b]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + 
					Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[1+3][i])/(3.0*cw) + Z_N[2][b]*conj(Yd[1])*conj(Z_D[1][i]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][i]*Z_N[2][a])/(32.0*pow(PI, 2));

				}
			}
		}
	}
	


    //mixed ->


	scalar_t Cp1_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 

				Cp1_mixed += -D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[1+3][i]*Z_D[1+3][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[5][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[2])*conj(Z_D[2][i]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[5][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[2])*conj(Z_D[2][j])) + Z_D[5][i]*Z_D[5][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][i]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j])))/(96.0*pow(PI, 2)) - D2mix*Z_D[5][j]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*conj(Z_D[1+3][i])/(8.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->


}

//BS_2

C_mix_bs_2_SUSY::C_mix_bs_2_SUSY() : WilsonCoefficient("C_BS_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 3},                       //m_s 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(3050305, 3131, 0, 1)
    };
}

double C_mix_bs_2_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C2_chargedhiggs=0.;

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][1])*conj(V_CKM[j][1]); 	/* NM: added 1eration dependence, V_CKM[ie][2] -> V_CKM[ie][1] */
        
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c - 2*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C2_gluino=0.;

	double D2g,D0g;

	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie] and Z_D[5][ie] -> Z_D[1+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[2][j]*pow(g_3, 4.)*conj(Z_D[1+3][i])*conj(Z_D[1+3][j])/(288.*pow(PI, 2));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	scalar_t C2_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 1eration dependence, Yd[2] -> Yd[1] and V_CKM[Ke][2] -> V_CKM[Ke][1] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   						
						C2_chargino += 0.;
						
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C2_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					C2_neutralino += D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[1+3][i])/(3.0*cw) + Z_N[2][b]*conj(Yd[1])*conj(Z_D[1][i]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[1+3][j])/(3.0*cw) + 
					Z_N[2][b]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/
					(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][a])/(32.0*pow(PI, 2));


				}
			}
		}
	}
	


    //mixed ->


	scalar_t C2_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 
				
				C2_mixed += D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[2][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][i]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*conj(Z_D[2][j]) + 3.0*Z_D[2][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - 
					Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[1+3][i]))/(48.0*pow(PI, 2)) + D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][a])*conj(Z_D[1+3][i])*conj(Z_D[1+3][j])/(48.0*pow(PI, 2));


    		}
		}
	}

    //higgs PInguin ->


}

//BS_2_tilde

C_mix_bs_2_tilde_SUSY::C_mix_bs_2_tilde_SUSY() : WilsonCoefficient("CT_BS_2", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 3},                       //m_s 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(3050305, 3232, 0, 1)
    };
}

double C_mix_bs_2_tilde_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  

	scalar_t Cp2_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][1])*conj(V_CKM[j][1]); 	/* NM: added 1eration dependence, V_CKM[ie][2] -> V_CKM[ie][1] */

		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c-2.*D0h)/(128.*pow(PI,2.)*pow(M_W,4.));

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t Cp2_gluino=0.;

	double D2g,D0g;

	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie] and Z_D[5][ie] -> Z_D[1+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[5][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[1][i])*conj(Z_D[1][j])/(288.*pow(PI, 2.));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	
	scalar_t Cp2_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 1eration dependence, Yd[2] -> Yd[1] and V_CKM[Ke][2] -> V_CKM[Ke][1] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						
						Cp2_chargino += 0.;
						
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t Cp2_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					Cp2_neutralino += D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
					Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][b]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - 
					conj(Z_N[1][b])*cw)*conj(Z_D[1][j])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][j])*conj(Z_N[2][b]))/(32.0*pow(PI, 2));


				}
			}
		}
	}
	


    //mixed ->


	scalar_t Cp2_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 

				Cp2_mixed += D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[5][i]*pow(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + 
				conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]), 2)*conj(Z_D[5][j]) + 3.0*Z_D[5][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))*conj(Z_D[1][i]))/(48.0*pow(PI, 2)) + 
				D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*conj(Z_D[1][i])*conj(Z_D[1][j])/(48.0*pow(PI, 2));


    		}
		}
	}

    //higgs PInguin ->



}

//B3

C_mix_bs_3_SUSY::C_mix_bs_3_SUSY() : WilsonCoefficient("C_BS_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 3},                       //m_s 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(3050305, 7171, 0, 1)
    };
}

double C_mix_bs_3_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  

	scalar_t C3_chargedhiggs=0.;

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][1])*conj(V_CKM[j][1]); 	/* NM: added 1eration dependence, V_CKM[ie][2] -> V_CKM[ie][1] */

		C3_chargedhiggs += 0.;

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C3_gluino=0.;

	double D2g,D0g;

	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie] and Z_D[5][ie] -> Z_D[1+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		C3_gluino += -D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[2][j]*pow(g_3,4.)*conj(Z_D[1+3][i])*conj(Z_D[1+3][j])/(96.*pow(PI, 2));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	
	scalar_t C3_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 1eration dependence, Yd[2] -> Yd[1] and V_CKM[Ke][2] -> V_CKM[Ke][1] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {

						
						C3_chargino += -D0ch*M_ch[a]*M_ch[b]*pow(V_CKM[k][2], 2)*Z_m[1][a]*Z_m[1][b]*Z_U[k][i]*Z_U[k][j]*(-Q_e*Z_p[0][a]*conj(Z_U[k][i])*swi + Yu[k]*Z_p[1][a]*conj(Z_U[k+3][i]))*(-Q_e*Z_p[0][b]*conj(Z_U[k][j])*swi + Yu[k]*Z_p[1][b]*conj(Z_U[k+3][j]))*pow(conj(V_CKM[k][1]), 2)*pow(conj(Yd[1]), 2)/(32.0*pow(PI, 2));
						

					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C3_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					
					C3_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*((-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][j]*Z_N[2][b]) - (-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][b]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][j]*Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[1+3][i])/(3.0*cw) + Z_N[2][b]*conj(Yd[1])*conj(Z_D[1][i]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][i]*Z_N[2][a])/(32.0*pow(PI, 2));


				}
			}
		}
	}
	


    //mixed ->


	scalar_t C3_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 

				C3_mixed += D0mix*M_ch0[a]*M_g*Z_D[2][i]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][i]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*conj(Z_D[2][j])/(48.0*pow(PI, 2)) - D0mix*M_ch0[a]*M_g*Z_D[2][j]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[1+3][i])/(48.0*pow(PI, 2)) + 
				D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
				Yd[2]*Z_D[5][j]*Z_N[2][a])*conj(Z_D[1+3][i])*conj(Z_D[1+3][j])/(48.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->

}

//BS3_tilde

C_mix_bs_3_tilde_SUSY::C_mix_bs_3_tilde_SUSY() : WilsonCoefficient("CT_BS_3", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 3},                       //m_s 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(3050305, 7272, 0, 1)
    };
}

double C_mix_bs_3_tilde_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  

	scalar_t Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][1])*conj(V_CKM[j][1]); 	/* NM: added 1eration dependence, V_CKM[ie][2] -> V_CKM[ie][1] */
        
		Cp3_chargedhiggs += 0.;
	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie] and Z_D[5][ie] -> Z_D[1+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[5][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[1][i])*conj(Z_D[1][j])/(96.*pow(PI, 2));
		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	scalar_t Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added 1eration dependence, Yd[2] -> Yd[1] and V_CKM[Ke][2] -> V_CKM[Ke][1] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						
						Cp3_chargino += -D0ch*M_ch[a]*M_ch[b]*pow(V_CKM[k][2], 2)*pow(Yd[2], 2)*(-Q_e*Z_U[k][i]*conj(Z_p[0][b])*swi + Z_U[k+3][i]*conj(Yu[k])*conj(Z_p[1][b]))*(-Q_e*Z_U[k][j]*conj(Z_p[0][a])*swi + Z_U[k+3][j]*conj(Yu[k])*conj(Z_p[1][a]))*pow(conj(V_CKM[k][1]), 2)*conj(Z_m[1][a])*conj(Z_m[1][b])*conj(Z_U[k][i])*conj(Z_U[k][j])/(32.0*pow(PI, 2));
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					Cp3_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*(-(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[1][j])/(2.0*cw*sw) + 
					conj(Yd[1])*conj(Z_D[1+3][j])*conj(Z_N[2][b])) + (-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][b])/(3.0*cw) + Yd[2]*Z_D[2][j]*conj(Z_N[2][b]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
					conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a])))*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + 
					Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][b]))/(32.0*pow(PI, 2));

				}
			}
		}
	}
	


    //mixed ->


	scalar_t Cp3_mixed=0.;


	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 

				Cp3_mixed += D0mix*M_ch0[a]*M_g*Z_D[5][i]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + 
				conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]), 2)*conj(Z_D[5][j])/(48.0*pow(PI, 2)) - D0mix*M_ch0[a]*M_g*Z_D[5][j]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + 
				conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))*conj(Z_D[1][i])/(48.0*pow(PI, 2)) + D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*conj(Z_D[1][i])*conj(Z_D[1][j])/(48.0*pow(PI, 2));
    		}
		}
	}

    //higgs PInguin ->


}

//B4

C_mix_bs_4_SUSY::C_mix_bs_4_SUSY() : WilsonCoefficient("C_BS_4", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 3},                       //m_s 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(3050305, 3132, 0, 1)
    };
}

double C_mix_bs_4_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  
	scalar_t C4_chargedhiggs=0.;

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][1])*conj(V_CKM[j][1]); 	/* NM: added 1eration dependence, V_CKM[ie][2] -> V_CKM[ie][1] */
        
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[i],2)*pow(m_u[j],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(PI,2.));

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C4_gluino=0.;

	double D2g,D0g;

	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie] and Z_D[5][ie] -> Z_D[1+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[1][i])*conj(Z_D[1+3][j])/(48.*pow(PI, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[1][i])*conj(Z_D[1+3][j]) + 11.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[1][j])*conj(Z_D[1+3][i]))/(72.*pow(PI, 2));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */

	scalar_t C4_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 1eration dependence, Yd[2] -> Yd[1] and V_CKM[Ke][2] -> V_CKM[Ke][1] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						
						C4_chargino += D2ch*pow(V_CKM[k][2], 2)*Yd[2]*Z_m[1][a]*Z_U[k][j]*(-Q_e*Z_p[0][b]*conj(Z_U[k][j])*swi + Yu[k]*Z_p[1][b]*conj(Z_U[k+3][j]))*(-Q_e*Z_U[k][i]*conj(Z_p[0][b])*swi + Z_U[k+3][i]*conj(Yu[k])*conj(Z_p[1][b]))*pow(conj(V_CKM[k][1]), 2)*conj(Yd[1])*conj(Z_m[1][a])*conj(Z_U[k][i])/(8.0*pow(PI, 2));

					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C4_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][j]*Z_N[2][b]) + (-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][b]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][j]*Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + 
					conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][b]))/(8.0*pow(PI, 2));


				}
			}
		}
	}
	


    //mixed ->


	scalar_t C4_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 

				C4_mixed += D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[2][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))*conj(Z_D[1+3][i]) + Z_D[5][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[1][i]))/(16.0*pow(PI, 2)) - 
				D2mix*pow(g_3, 2)*(-3.0*Z_D[2][j]*Z_D[5][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][i]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a])) - 3.0*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[5][i])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[2])*conj(Z_D[2][i]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][a])*conj(Z_D[1][j])*conj(Z_D[1+3][i]))/(24.0*pow(PI, 2)) - 
				D2mix*pow(g_3, 2)*(-Z_D[2][j]*Z_D[5][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a])) - (-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[1][j])*conj(Z_D[1+3][i]))/(24.0*pow(PI, 2)) - 
				D2mix*pow(g_3, 2)*(Z_D[2][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*conj(Z_D[1][i]) + Z_D[5][j]*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + 
				Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))*conj(Z_D[1+3][i]))/(24.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->




	scalar_t delta_d[6][6],delta_d_LL[3][3],delta_d_LR[3][3],delta_d_RL[3][3],delta_d_RR[3][3];
	scalar_t delta_u[6][6],delta_u_LL[3][3],delta_u_LR[3][3],delta_u_RL[3][3],delta_u_RR[3][3];
	
	double m_av = (M_U[0]+M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_D[0]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double M_A=src.at({ParameterType::BSM, "MASS", 36})->get_val();
	double A_t=src.at({ParameterType::BSM, "AU", LhaID(3,3)})->get_val(); //A_t
	double mu = src.at({ParameterType::BSM, "HMIX", 1})->get_val();
	scalar_t M_2 = src.at({ParameterType::BSM, "MSOFT", 2})->get_val();
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	scalar_t x_2 = pow(abs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*PI) ;
	double alpha_2 = pow(g_2,2.)/(4.*PI);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*PI*pow(m_av,2));
	
	
	scalar_t V_tb = V_CKM[2][2]; 
	scalar_t V_tq = V_CKM[2][1];
	scalar_t a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*PI*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	scalar_t a2 = -alpha_s*pow(M_g,2)*delta_d_LL[2][1]*delta_d_RR[2][1]*pow(h1(x_g),2);
	scalar_t a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[2][1]*V_tb*conj(V_tq)/(pow(M_W,2));
	scalar_t a4 = alpha_2*M_2*M_g*delta_u_LL[2][1]*delta_d_RR[2][1]*h1(x_g)*h4(x_2,x_g);

	scalar_t C4_higgspenguin = a1*(a2+a3+a4);

}

//BS_5

C_mix_bs_5_SUSY::C_mix_bs_5_SUSY() : WilsonCoefficient("C_BS_5", GroupMapper::str(WGroup::MESON_MIXING) + "_MATCH") {
    matching_info[QCDOrder::LO] = {
        {
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 4},           //mass_c_muW_mcrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", LhaID(5, 1)}, //mass_b_muW_mbrun
            {ParameterType::WILSON, "WPARAM_MATCH_SM", 6},           //mass_t_muW_mbrun
            {ParameterType::SM, "MASS", 24},  
            {ParameterType::BSM, "MASS", 37},                     // M_H
            {ParameterType::SM, "MASS", 3},                       //m_s 
            {ParameterType::SM, "MASS", 2},                       //m_u
			{ParameterType::BSM, "MASS", 1000001},
			{ParameterType::BSM, "MASS", 1000002},
			{ParameterType::BSM, "MASS", 1000003},
			{ParameterType::BSM, "MASS", 1000004},
			{ParameterType::BSM, "MASS", 1000005},
			{ParameterType::BSM, "MASS", 1000006},
			{ParameterType::BSM, "MASS", 2000001},
			{ParameterType::BSM, "MASS", 2000002},
			{ParameterType::BSM, "MASS", 2000003},
			{ParameterType::BSM, "MASS", 2000004},
			{ParameterType::BSM, "MASS", 2000005},
			{ParameterType::BSM, "MASS", 2000006},
			{ParameterType::BSM, "MASS", 1000024},
			{ParameterType::BSM, "MASS", 1000037},
			{ParameterType::BSM, "MASS", 1000022},
			{ParameterType::BSM, "MASS", 1000023},
			{ParameterType::BSM, "MASS", 1000025},
			{ParameterType::BSM, "MASS", 1000035},
			{ParameterType::BSM, "MASS", 36},
			{ParameterType::BSM, "MSOFT", 2},
			{ParameterType::BSM, "HMIX", 1},
			{ParameterType::BSM, "AU", LhaID(3,3)},
			{ParameterType::SM, "GAUGE", 1},						//gp 
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
			{ParameterType::SM, "EW_SCALE", 1}
        },
        compute_LO,
        LhaID(3050305, 7172, 0, 1)
    };
}

double C_mix_bs_5_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

	double mu_W = src.at({ParameterType::WILSON, "EW_SCALE", 37})->get_val();
    double M_H=src.at({ParameterType::BSM, "MASS", 37})->get_val();
	double M_W=src.at({ParameterType::SM, "MASS", 24})->get_val();
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	std::array<std::array<scalar_t, 3>, 3> V_CKM {};
	double m_q = src.at({ParameterType::SM, "MASS", 3})->get_val();
	
	double m_b= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", {5, 1}})->get_val();
	double g_2=src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
	double m_u[4],m_u_pow_2[4];

    for (int i = 0; i<3; ++i) {
        for (int j = 0; j<3; j++) {
            V_CKM[i][j] = src.at({ParameterType::SM, "VCKM", LhaID(i, j)})->get_val();
        }
    }
	m_u[1]= src.at({ParameterType::SM, "MASS", 2})->get_val();
    m_u[2]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 4})->get_val();
	m_u[3]= src.at({ParameterType::WILSON, "WPARAM_MATCH_SM", 6})->get_val();


	for(int i =0;i<3;++i) {
        m_u_pow_2[i]=pow(m_u[i],2.);
    }
  

	scalar_t C5_chargedhiggs=0.;

	
	double D0h,D0h_c,D2h,D2h_c;
	scalar_t CKM_product;


    for(int i = 0; i<3; i++) for(int j=1; j<3; j++) {
		D0h = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[i],m_u_pow_2[j],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[i][3]*V_CKM[j][3]*conj(V_CKM[i][1])*conj(V_CKM[j][1]); 	/* NM: added 1eration dependence, V_CKM[ie][2] -> V_CKM[ie][1] */
        
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[i],2.)*(D2h_c-2.*D2h)/(32.*pow(PI,2.)*pow(M_W,4.));

	}

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	// double g_3=sqrt(4.*PI*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_3= sqrt(4.*PI*QCDHelper::alpha_s(mu_W)); //TODO : check pole or running
	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= src.at({ParameterType::BSM, "DSQMIX", LhaID(j+1, i+1)})->get_val(); //TODO: deal with this group
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	scalar_t C5_gluino=0.;

	double D2g,D0g;

	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie] and Z_D[5][ie] -> Z_D[1+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);

		C5_gluino += -D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[1][i])*conj(Z_D[1+3][j])/(144.*pow(PI, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[1][i])*conj(Z_D[1+3][j]) + 3.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[1][j])*conj(Z_D[1+3][i]))/(72.*pow(PI, 2));

		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(src.at({ParameterType::SM, "GAUGE", 2})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));
 	double Q_e = (src.at({ParameterType::SM, "GAUGE", 2})->get_val())*sw;
	double swi=1./sw;

	M_ch[0]=src.at({ParameterType::BSM, "MASS", 1000024})->get_val();
	M_ch[1]=src.at({ParameterType::BSM, "MASS", 1000037})->get_val();

	M_U[0]=src.at({ParameterType::BSM, "MASS", 1000002})->get_val();
	M_U[1]=src.at({ParameterType::BSM, "MASS", 1000004})->get_val();
	M_U[2]=src.at({ParameterType::BSM, "MASS", 1000006})->get_val();
	M_U[3]=src.at({ParameterType::BSM, "MASS", 2000002})->get_val();
	M_U[4]=src.at({ParameterType::BSM, "MASS", 2000004})->get_val();
	M_U[5]=src.at({ParameterType::BSM, "MASS", 2000006})->get_val();

	for(int i=0; i<2; ++i){
		M_ch_pow_2[i]=pow(M_ch[i],2);
	}
	for(int i=0; i<6; ++i){
		M_U_pow_2[i]=pow(M_U[i],2);
	}
	for(int i=0; i <2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_p[i][j]=conj(src.at({ParameterType::BSM, "VMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=2; i<2; ++i) {
		for(int j=0; j<2; ++j) {
			Z_m[i][j]=conj(src.at({ParameterType::BSM, "UMIX", LhaID(j+1, i+1)})->get_val());
		}
	} /* NM: conversion from SLHA2 convention */
	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			Z_U[i][j]=conj(src.at({ParameterType::BSM, "USQMIX", LhaID(j+1, i+1)})->get_val()); //TODO: deal with this group
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double m_b=QCDHelper::msbar_mass(5, mu_W, MassType::MSBAR); /* NM: running mass */
	double m_t=QCDHelper::msbar_mass(6, mu_W, MassType::MSBAR); /* NM: running mass */

	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*QCDHelper::msbar_mass(4, mu_W, MassType::POLE); //TODO : check this to be sure

	Yu[2] = common*m_t; /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*m_b; /* NM: running mass */
	

	scalar_t C5_chargino=0.;

  
	double D0ch,D2ch;

	/* NM: added 1eration dependence, Yd[2] -> Yd[1] and V_CKM[Ke][2] -> V_CKM[Ke][1] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {

						C5_chargino += -D0ch*M_ch[a]*M_ch[b]*pow(V_CKM[k][2], 2)*Yd[2]*Z_m[1][b]*Z_U[k][i]*(-Q_e*Z_p[0][b]*conj(Z_U[k][j])*swi + Yu[k]*Z_p[1][b]*conj(Z_U[k+3][j]))*(-Q_e*Z_U[k][j]*conj(Z_p[0][a])*swi + Z_U[k+3][j]*conj(Yu[k])*conj(Z_p[1][a]))*pow(conj(V_CKM[k][1]), 2)*conj(Yd[1])*conj(Z_m[1][a])*conj(Z_U[k][i])/(16.0*pow(PI, 2));
					
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[4],M_ch0_pow_2[4],M_D[6],M_D_pow_2[6];
	scalar_t Z_N[4][4];

	double cw=cos(atan(src.at({ParameterType::SM, "GAUGE", 1})->get_val()/src.at({ParameterType::SM, "GAUGE", 2})->get_val()));

 	double otherc = sqrt(2.)/v1;

	std::array<double,4> temp_ch0 = {src.at({ParameterType::BSM, "MASS", 1000022})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000023})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000025})->get_val(),
					 src.at({ParameterType::BSM, "MASS", 1000035})->get_val()};

	M_ch0[0]=fabs(temp_ch0[0]);
	M_ch0[1]=fabs(temp_ch0[1]);
	M_ch0[2]=fabs(temp_ch0[2]);
	M_ch0[3]=fabs(temp_ch0[3]);

	for(int i=0; i<4; ++i){
		M_ch0_pow_2[i]=pow(M_ch0[i],2);
	}
	
	
	for(int i=0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}


	
	for(int i=0; i<4; ++i) {
		for(int j=0; j<4; ++j) {
			Z_N[i][j] = conj(src.at({ParameterType::BSM, "NMIX", LhaID(i+1, j+1)})->get_val()); //TODO : i,j or j,i like the others ?
		}
	}
	

	for(int i=0; i<4; ++i){
		if(temp_ch0[i]<0.) {
			for(int j=0; j<4; ++j) {
				Z_N[i][j]*=I;
			}
		}
	} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C5_neutralino=0.;

	double D0ne,D2ne;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				for (int b=0; b<4; ++b) {
					D0ne = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);
					D2ne = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_ch0_pow_2[b]);

					C5_neutralino += -D0ne*M_ch0[a]*M_ch0[b]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[1+3][i])/(3.0*cw) + 
					Z_N[2][b]*conj(Yd[1])*conj(Z_D[1][i]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][b]*sw/3.0 - Z_N[1][b]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][b])*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
					conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))/(16.0*pow(PI, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
					Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][b]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][b]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0- Z_N[1][a]*cw)/(2.0*cw*sw) + 
					Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*(conj(Z_N[0][b])*sw/3.0 - conj(Z_N[1][b])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][b]))/(8.0*pow(PI, 2));

				}
			}
		}
	}
	


    //mixed ->

	scalar_t C5_mixed=0.;



	double M_g=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double M_g_pow_2 = pow(M_g,2.);

	// } /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added 1eration dependence, Z_D[2][ie] -> Z_D[1][ie], Z_D[5][ie] -> Z_D[1+3][ie] and Yd[2] -> Yd[1] */

	for(int i=0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a=0; a<4; ++a) {
				D0mix = D0(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2);
				D2mix = D2p(M_D_pow_2[i],M_D_pow_2[j],M_ch0_pow_2[a],M_g_pow_2); 

				C5_mixed += -D0mix*M_ch0[a]*M_g*pow(g_3, 2)*(Z_D[2][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))*conj(Z_D[1+3][i]) + Z_D[5][j]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[1][i]))/(48.0*pow(PI, 2)) + 
				D2mix*pow(g_3, 2)*(-Z_D[2][j]*Z_D[5][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][i])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][i]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a])) - (-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[5][i])/(3*cw) + 
				Z_N[2][a]*conj(Yd[2])*conj(Z_D[2][i]))*(-sqrt(2)*Q_e*Z_D[2][j]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][j]*Z_N[2][a])*conj(Z_D[1][j])*conj(Z_D[1+3][i]))/(24.0*pow(PI, 2)) + 
				D2mix*pow(g_3, 2)*(-Z_D[2][j]*Z_D[5][i]*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a])) - (-sqrt(2)*Q_e*Z_D[5][j]*conj(Z_N[0][a])/(3.0*cw) + 
				Yd[2]*Z_D[2][j]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*conj(Z_D[1][j])*conj(Z_D[1+3][i]))/(8.0*pow(PI, 2)) + 
				D2mix*pow(g_3, 2)*(Z_D[2][j]*(-sqrt(2)*Q_e*Z_D[5][i]*conj(Z_N[0][a])/(3.0*cw) + Yd[2]*Z_D[2][i]*conj(Z_N[2][a]))*(-sqrt(2)*Q_e*Z_N[0][a]*conj(Z_D[1+3][j])/(3.0*cw) + 
				Z_N[2][a]*conj(Yd[1])*conj(Z_D[1][j]))*conj(Z_D[1][i]) + Z_D[5][j]*(-sqrt(2)*Q_e*Z_D[2][i]*(Z_N[0][a]*sw/3.0 - Z_N[1][a]*cw)/(2.0*cw*sw) + Yd[2]*Z_D[5][i]*Z_N[2][a])*(-sqrt(2)*Q_e*(conj(Z_N[0][a])*sw/3.0 - 
				conj(Z_N[1][a])*cw)*conj(Z_D[1][i])/(2.0*cw*sw) + conj(Yd[1])*conj(Z_D[1+3][i])*conj(Z_N[2][a]))*conj(Z_D[1+3][i]))/(8.0*pow(PI, 2));

    		}
		}
	}

    //higgs PInguin ->


}