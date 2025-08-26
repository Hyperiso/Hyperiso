#include "MesonMixingWilson_SUSY.h"

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
            {ParameterType::SM, "MASS", 3},                       //m_s
			{ParameterType::BSM, "MASS", 1000001}
			{ParameterType::BSM, "MASS", 1000002}
			{ParameterType::BSM, "MASS", 1000003}
			{ParameterType::BSM, "MASS", 1000004}
			{ParameterType::BSM, "MASS", 1000005}
			{ParameterType::BSM, "MASS", 1000006}
			{ParameterType::BSM, "MASS", 2000001}
			{ParameterType::BSM, "MASS", 2000002}
			{ParameterType::BSM, "MASS", 2000003}
			{ParameterType::BSM, "MASS", 2000004}
			{ParameterType::BSM, "MASS", 2000005}
			{ParameterType::BSM, "MASS", 2000006}
			{ParameterType::BSM, "MASS", 1000024}
			{ParameterType::BSM, "MASS", 1000037}
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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta = src.at({ParameterType::BSM, "HMIX", 2})->get_val();
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

    //gluino ->


	double M_D[6],M_D_pow_2[6],dm[6]; 
	scalar_t Z_D[6][6];
	double Mg=src.at({ParameterType::BSM, "MASS", 1000021})->get_val();
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[0]=src.at({ParameterType::BSM, "MASS", 1000001})->get_val();
	M_D[1]=src.at({ParameterType::BSM, "MASS", 1000003})->get_val();
	M_D[2]=src.at({ParameterType::BSM, "MASS", 1000005})->get_val();
	M_D[3]=src.at({ParameterType::BSM, "MASS", 2000001})->get_val();
	M_D[4]=src.at({ParameterType::BSM, "MASS", 2000003})->get_val();
	M_D[5]=src.at({ParameterType::BSM, "MASS", 2000005})->get_val();

	for(int i = 0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
			Z_D[i][j]= param->sD_mix[j][i];
		}
	} // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(int i = 0; i<6; ++i) {
		M_D_pow_2[i]=pow(M_D[i],2);
	}

	scalar_t C1_gluino=0.;
	scalar_t C2_gluino=0.;
	scalar_t C3_gluino=0.;
	scalar_t C4_gluino=0.;
	scalar_t C5_gluino=0.;
	scalar_t Cp1_gluino=0.;
	scalar_t Cp2_gluino=0.;
	scalar_t Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(int i =0; i<6; ++i) {
		for(int j = 0; j<6; ++j) {
		D2g = D2p(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[i],M_D_pow_2[j],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[2][i]*Z_D[2][j]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][i])*conj(Z_D[gen][j])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[2][j]*pow(g_3, 4.)*conj(Z_D[gen+3][i])*conj(Z_D[gen+3][j])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[2][j]*pow(g_3,4.)*conj(Z_D[gen+3][i])*conj(Z_D[gen+3][j])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[gen][i])*conj(Z_D[gen+3][j])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[gen][i])*conj(Z_D[gen+3][j]) + 11.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[gen][j])*conj(Z_D[gen+3][i]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[2][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[gen][i])*conj(Z_D[gen+3][j])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[gen][i])*conj(Z_D[gen+3][j]) + 3.*Z_D[2][i]*Z_D[5][j]*conj(Z_D[gen][j])*conj(Z_D[gen+3][i]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[5][i]*Z_D[5][j]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][i])*conj(Z_D[gen+3][j])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[5][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[gen][i])*conj(Z_D[gen][j])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[5][i]*Z_D[5][j]*pow(g_3,4.)*conj(Z_D[gen][i])*conj(Z_D[gen][j])/(96.*pow(pi, 2));
		}
	}


    //chargino ->

	
	double M_ch[2],M_ch_pow_2[2],M_U[6],M_U_pow_2[6];
	scalar_t Z_p[2][2],Z_m[2][2],Z_U[6][6],V_CKM[3][3];
	scalar_t Yd[3],Yu[3];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
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
			Z_U[i][j]=conj(param->sU_mix[j][i]);
		}
	} /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(src.at({ParameterType::BSM, "EXTPAR", 25})->get_val());
	v1 = 2.*(src.at({ParameterType::SM, "MASS", 24})->get_val())*cos(beta)/src.at({ParameterType::SM, "GAUGE", 2})->get_val();
	v2 = v1*tan(beta);
  
	double mc = src.at({ParameterType::SM, "MASS", 4})->get_val();
	
	double common = sqrt(2.)/v2;
	Yu[0] = common*src.at({ParameterType::SM, "MASS", 2})->get_val();
	Yu[1] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[2] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[0] = otherc*src.at({ParameterType::SM, "MASS", 1})->get_val();
	Yd[1] = otherc*src.at({ParameterType::SM, "MASS", 3})->get_val();
	Yd[2] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(int i = 0; i<3; ++i) {
		for(int j = 0; j<3; ++j) {
			V_CKM[i][j]=src.at({ParameterType::SM, "VCKM", LhaID(i,j)})->get_val();
		}
	}
	
	scalar_t C1_chargino=0.;
	scalar_t C2_chargino=0.;
	scalar_t C3_chargino=0.;
	scalar_t C4_chargino=0.;
	scalar_t C5_chargino=0.;
	scalar_t Cp1_chargino=0.;
	scalar_t Cp2_chargino=0.;
	scalar_t Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(int i = 0; i<6; ++i) {
		for(int j=0; j<6; ++j) {
			for(int a =0; a<2; ++a) {
				for (int b=0; b<2; ++b) {
					D0ch = D0(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]);
					D2ch = D2p(M_U_pow_2[i],M_U_pow_2[j],M_ch_pow_2[a],M_ch_pow_2[b]); 
					for(int k=0; k<3; ++k) {
						   
						
						C1_chargino += -D2ch*pow(V_CKM[k][2], 2)*(-Q_e*Z_p[0][a]*conj(Z_U[k][i])*swi + Yu[k]*Z_p[1][a]*conj(Z_U[k+3][i]))*(-Q_e*Z_p[0][b]*conj(Z_U[k][j])*swi + Yu[k]*Z_p[1][b]*conj(Z_U[k+3][j]))*(-Q_e*Z_U[k][i]*conj(Z_p[0][b])*swi + Z_U[k+3][i]*conj(Yu[k])*conj(Z_p[1][b]))*(-Q_e*Z_U[k][j]*conj(Z_p[0][a])*swi + Z_U[k+3][j]*conj(Yu[k])*conj(Z_p[1][a]))*pow(conj(V_CKM[k][gen]), 2)/(32.0*pow(pi, 2));
						
						C2_chargino += 0.;
						
						C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
						
						C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
						
						C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
					
						Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
						
						Cp2_chargino += 0.;
						
						Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
					}
				}
			}
		}
	}
	


    //neutralino ->


	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	scalar_t Z_D[7][7],Z_N[5][5];
	scalar_t Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	// M_D[1]=param->mass_dnl;
	// M_D[2]=param->mass_stl;
	// M_D[3]=param->mass_b1;
	// M_D[4]=param->mass_dnr;
	// M_D[5]=param->mass_str;
	// M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	scalar_t C1_neutralino=0.;
	scalar_t C2_neutralino=0.;
	scalar_t C3_neutralino=0.;
	scalar_t C4_neutralino=0.;
	scalar_t C5_neutralino=0.;
	scalar_t Cp1_neutralino=0.;
	scalar_t Cp2_neutralino=0.;
	scalar_t Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	


    //mixed ->


	scalar_t C1_mixed=0.;
	scalar_t C2_mixed=0.;
	scalar_t C3_mixed=0.;
	scalar_t C4_mixed=0.;
	scalar_t C5_mixed=0.;
	scalar_t Cp1_mixed=0.;
	scalar_t Cp2_mixed=0.;
	scalar_t Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	// M_D[1]=param->mass_dnl;
	// M_D[2]=param->mass_stl;
	// M_D[3]=param->mass_b1;
	// M_D[4]=param->mass_dnr;
	// M_D[5]=param->mass_str;
	// M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }

    //higgs pinguin ->


	// double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	// M_D[1]=param->mass_dnl;
	// M_D[2]=param->mass_stl;
	// M_D[3]=param->mass_b1;
	// M_D[4]=param->mass_dnr;
	// M_D[5]=param->mass_str;
	// M_D[6]=param->mass_b2;
	// M_U[1]=param->mass_upl;
	// M_U[2]=param->mass_chl;
	// M_U[3]=param->mass_t1;
	// M_U[4]=param->mass_upr;
	// M_U[5]=param->mass_chr;
	// M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	scalar_t delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	scalar_t delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	scalar_t V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	scalar_t V_tb = V_CKM[3][3]; 
	scalar_t V_tq = V_CKM[3][gen];
	scalar_t a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	scalar_t a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	scalar_t a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	scalar_t a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	scalar_t C4_higgspenguin = a1*(a2+a3+a4);

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}

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

double C_mix_bd_1_SUSY::compute_LO(const std::unordered_map<ParamId, std::shared_ptr<Parameter>>& src) {

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
	double tbeta=param->tan_beta;
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

    //gluino ->
    
    int ie,je;
	for(ie=1;ie<=5;ie++) CM_gluino[ie];
	for(ie=1;ie<=3;ie++) CMp_gluino[ie]=0.;

	double M_D[7],M_D_pow_2[7],dm[7]; 
	double complex Z_D[7][7];
	double Mg=param->mass_gluino;
	double Mg_pow_2 = pow(Mg,2);
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; // inverse matrix, because in SLHA2 the second index denotes quark flavour (dl,sl,bl,dr,sr,br)
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);

	double complex C1_gluino=0.;
	double complex C2_gluino=0.;
	double complex C3_gluino=0.;
	double complex C4_gluino=0.;
	double complex C5_gluino=0.;
	double complex Cp1_gluino=0.;
	double complex Cp2_gluino=0.;
	double complex Cp3_gluino=0.;
	double D2g,D0g;

	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie] and Z_D[5][ie] -> Z_D[gen+3][ie] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		D2g = D2p(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		D0g = D0(M_D_pow_2[ie],M_D_pow_2[je],Mg_pow_2,Mg_pow_2);
		C1_gluino += -Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*(D0g*pow(Mg, 2) + 11.*D2g)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(144.*pow(pi, 2));
		C2_gluino += -17.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3, 4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(288.*pow(pi, 2));
		C3_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[3][je]*pow(g_3,4.)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(96.*pow(pi, 2));
		C4_gluino += -7.*D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(48.*pow(pi, 2)) + D2g*pow(g_3,4.)*(6.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 11.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		C5_gluino += -D0g*pow(Mg, 2)*Z_D[3][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2)) + 5.*D2g*pow(g_3,4.)*(-2.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][ie])*conj(Z_D[gen+3][je]) + 3.*Z_D[3][ie]*Z_D[6][je]*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(72.*pow(pi, 2));
		Cp1_gluino += -Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*(D0g*pow(Mg, 2.) + 11.*D2g)*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(144.*pow(pi, 2.));
		Cp2_gluino += -17.*D0g*pow(Mg, 2.)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(288.*pow(pi, 2.));
		Cp3_gluino += -D0g*pow(Mg, 2)*Z_D[6][ie]*Z_D[6][je]*pow(g_3,4.)*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(96.*pow(pi, 2));
	}
	CM_gluino[1]=C1_gluino;
	CM_gluino[2]=C2_gluino;
	CM_gluino[3]=C3_gluino;
	CM_gluino[4]=C4_gluino;
	CM_gluino[5]=C5_gluino;
	CMp_gluino[1]=Cp1_gluino;
	CMp_gluino[2]=Cp2_gluino;
	CMp_gluino[3]=Cp3_gluino;


    //chargino ->

    int ie,je,ae,be,Ke;
	for(ie=1;ie<=5;ie++) CM_chargino[ie];
	for(ie=1;ie<=3;ie++) CMp_chargino[ie]=0.;
	
	double M_ch[3],M_ch_pow_2[3],M_U[7],M_U_pow_2[7];
	double complex Z_p[3][3],Z_m[3][3],Z_U[7][7],V_CKM[4][4];
	double complex Yd[4],Yu[4];
 	double sw=sin(atan(param->gp/param->g2));
 	double Q_e = (param->g2)*sw;
	double swi=1./sw;

	M_ch[1]=param->mass_cha1;
	M_ch[2]=param->mass_cha2;

	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	for(ie=1;ie<=2;ie++){M_ch_pow_2[ie]=pow(M_ch[ie],2);}
	for(ie=1;ie<=6;ie++){M_U_pow_2[ie]=pow(M_U[ie],2);}
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_p[ie][je]=conj(param->charg_Vmix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) Z_m[ie][je]=conj(param->charg_Umix[je][ie]); /* NM: conversion from SLHA2 convention */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]=conj(param->sU_mix[je][ie]); /* NM: conversion from SLHA2 convention */

	double v1,v2,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
	v2 = v1*tan(beta);
  
	double common = sqrt(2.)/v2;
	Yu[1] = common*param->mass_u;
	Yu[2] = common*running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	Yu[3] = common*running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex C1_chargino=0.;
	double complex C2_chargino=0.;
	double complex C3_chargino=0.;
	double complex C4_chargino=0.;
	double complex C5_chargino=0.;
	double complex Cp1_chargino=0.;
	double complex Cp2_chargino=0.;
	double complex Cp3_chargino=0.;
  
	double D0ch,D2ch;

	/* NM: added generation dependence, Yd[2] -> Yd[gen] and V_CKM[Ke][2] -> V_CKM[Ke][gen] */
  
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=2;ae++) for (be=1;be<=2;be++) for(Ke=1;Ke<=3;Ke++)
	{
		D0ch = D0(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);
		D2ch = D2p(M_U_pow_2[ie],M_U_pow_2[je],M_ch_pow_2[ae],M_ch_pow_2[be]);    
		
		C1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)/(32.0*pow(pi, 2));
		
		C2_chargino += 0.;
		
		C3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*(-Q_e*Z_p[1][ae]*conj(Z_U[Ke][ie])*swi + Yu[Ke]*Z_p[2][ae]*conj(Z_U[Ke+3][ie]))*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)/(32.0*pow(pi, 2));
		
		C4_chargino += D2ch*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][ae]*Z_U[Ke][je]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(8.0*pow(pi, 2));
		
		C5_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*Yd[3]*Z_m[2][be]*Z_U[Ke][ie]*(-Q_e*Z_p[1][be]*conj(Z_U[Ke][je])*swi + Yu[Ke]*Z_p[2][be]*conj(Z_U[Ke+3][je]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Yd[gen])*conj(Z_m[2][ae])*conj(Z_U[Ke][ie])/(16.0*pow(pi, 2));
    
		Cp1_chargino += -D2ch*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*Z_m[2][ae]*Z_m[2][be]*Z_U[Ke][ie]*Z_U[Ke][je]*pow(conj(V_CKM[Ke][gen]), 2)*pow(conj(Yd[gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
		
		Cp2_chargino += 0.;
		
		Cp3_chargino += -D0ch*M_ch[ae]*M_ch[be]*pow(V_CKM[Ke][3], 2)*pow(Yd[3], 2)*(-Q_e*Z_U[Ke][ie]*conj(Z_p[1][be])*swi + Z_U[Ke+3][ie]*conj(Yu[Ke])*conj(Z_p[2][be]))*(-Q_e*Z_U[Ke][je]*conj(Z_p[1][ae])*swi + Z_U[Ke+3][je]*conj(Yu[Ke])*conj(Z_p[2][ae]))*pow(conj(V_CKM[Ke][gen]), 2)*conj(Z_m[2][ae])*conj(Z_m[2][be])*conj(Z_U[Ke][ie])*conj(Z_U[Ke][je])/(32.0*pow(pi, 2));
	}
	
	CM_chargino[1]=C1_chargino;
	CM_chargino[2]=C2_chargino;
	CM_chargino[3]=C3_chargino;
	CM_chargino[4]=C4_chargino;
	CM_chargino[5]=C5_chargino;
	CMp_chargino[1]=Cp1_chargino;
	CMp_chargino[2]=Cp2_chargino;
	CMp_chargino[3]=Cp3_chargino;


    //neutralino ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_neutralino[ie];
	for(ie=1;ie<=3;ie++) CMp_neutralino[ie]=0.;

	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];
	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	
	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je] = param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */
	
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je] = param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double complex C1_neutralino=0.;
	double complex C2_neutralino=0.;
	double complex C3_neutralino=0.;
	double complex C4_neutralino=0.;
	double complex C5_neutralino=0.;
	double complex Cp1_neutralino=0.;
	double complex Cp2_neutralino=0.;
	double complex Cp3_neutralino=0.;
	double D0ne,D2ne;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++) for (be=1;be<=4;be++)
	{
		D0ne = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		D2ne = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_ch0_pow_2[be]);
		C1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		C2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) - (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		C4_neutralino += D2ne*((-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be]) + (-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));
		C5_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][be]*sw/3.0 - Z_N[2][be]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][be])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))/(16.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0- Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(8.0*pow(pi, 2));

		Cp1_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))/(64.0*pow(pi, 2)) - D2ne*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_N[1][be]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][be]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])/(32.0*pow(pi, 2));
		Cp2_neutralino += D0ne*M_ch0[ae]*M_ch0[be]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
		Cp3_neutralino += -D0ne*M_ch0[ae]*M_ch0[be]*(-(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][je])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][je])*conj(Z_N[3][be])) + (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][be])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][be]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])))*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][be])*sw/3.0 - conj(Z_N[2][be])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][be]))/(32.0*pow(pi, 2));
	}
	
	CM_neutralino[1]=C1_neutralino;
	CM_neutralino[2]=C2_neutralino;
	CM_neutralino[3]=C3_neutralino;
	CM_neutralino[4]=C4_neutralino;
	CM_neutralino[5]=C5_neutralino;
	CMp_neutralino[1]=Cp1_neutralino;
	CMp_neutralino[2]=Cp2_neutralino;
	CMp_neutralino[3]=Cp3_neutralino;


    //mixed ->

    int ie,je,ae,be;
	for(ie=1;ie<=5;ie++) CM_mixed[ie];
	for(ie=1;ie<=3;ie++) CMp_mixed[ie]=0.;

	double complex C1_mixed=0.;
	double complex C2_mixed=0.;
	double complex C3_mixed=0.;
	double complex C4_mixed=0.;
	double complex C5_mixed=0.;
	double complex Cp1_mixed=0.;
	double complex Cp2_mixed=0.;
	double complex Cp3_mixed=0.;

	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */

	double sw=sin(atan(param->gp/param->g2));
	double cw=cos(atan(param->gp/param->g2));
	double Q_e = param->g2*sw;
	double M_g=param->mass_gluino;
	double M_g_pow_2 = pow(M_g,2.);
	double M_ch0[5],M_ch0_pow_2[5],M_D[7],M_D_pow_2[7];
	double complex Z_D[7][7],Z_N[5][5];
	double complex Yd[4];

	double v1,beta;
	beta = atan(param->tan_beta);
	v1 = 2.*(param->mass_W)*cos(beta)/param->g2;
 	double otherc = sqrt(2.)/v1;
	Yd[1] = otherc*param->mass_d;
	Yd[2] = otherc*param->mass_s;
	Yd[3] = otherc*running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	M_ch0[1]=fabs(param->mass_neut[1]);
	M_ch0[2]=fabs(param->mass_neut[2]);
	M_ch0[3]=fabs(param->mass_neut[3]);
	M_ch0[4]=fabs(param->mass_neut[4]);
	for(ie=1;ie<=4;ie++){M_ch0_pow_2[ie]=pow(M_ch0[ie],2);}
	
	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;

	for(ie=1;ie<=6;ie++) M_D_pow_2[ie]=pow(M_D[ie],2);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie]; /* NM: conversion from SLHA2 convention */

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=4;ie++) for(je=1;je<=4;je++) Z_N[ie][je]=param->neut_mix[ie][je];
	for(ie=1;ie<=4;ie++){if(param->mass_neut[ie]<0.) for(je=1;je<=4;je++) Z_N[ie][je]*=I;} /* NM: in Buras, neutralino masses defined positive, so changed the SLHA2 neutralino mixing matrix when negative neutralino mass */

	double D0mix,D2mix;
	
	/* NM: added generation dependence, Z_D[2][ie] -> Z_D[gen][ie], Z_D[5][ie] -> Z_D[gen+3][ie] and Yd[2] -> Yd[gen] */

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) for(ae=1;ae<=4;ae++)
	{
		D0mix = D0(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2);
		D2mix = D2p(M_D_pow_2[ie],M_D_pow_2[je],M_ch0_pow_2[ae],M_g_pow_2); 
		
		C1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen][ie]*Z_D[gen][je]*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][ie])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][ie])*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[3][je])/(2.0*cw*sw) + conj(Yd[3])*conj(Z_D[6][je])*conj(Z_N[3][ae])) + Z_D[3][ie]*Z_D[3][je]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2))/(96.0*pow(pi, 2)) - D2mix*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(8.0*pow(pi, 2));
		C2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je]) + 3.0*Z_D[3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[3][ie]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[3][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[3][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen+3][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen+3][ie])*conj(Z_D[gen+3][je])/(48.0*pow(pi, 2));
		C4_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(16.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-3.0*Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - 3.0*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) - D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2));
		C5_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_D[3][je]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][je]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(24.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(-Z_D[3][je]*Z_D[6][ie]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae])) - (-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*conj(Z_D[gen][je])*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2)) + D2mix*pow(g_3, 2)*(Z_D[3][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen][ie]) + Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[3][ie]*(Z_N[1][ae]*sw/3.0 - Z_N[2][ae]*cw)/(2.0*cw*sw) + Yd[3]*Z_D[6][ie]*Z_N[3][ae])*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen+3][ie]))/(8.0*pow(pi, 2));
		Cp1_mixed += -D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[gen+3][ie]*Z_D[gen+3][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[6][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[3])*conj(Z_D[3][je])) + Z_D[6][ie]*Z_D[6][je]*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][ie])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][ie]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je])))/(96.0*pow(pi, 2)) - D2mix*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_N[1][ae]*conj(Z_D[gen+3][je])/(3.0*cw) + Z_N[3][ae]*conj(Yd[gen])*conj(Z_D[gen][je]))*conj(Z_D[gen+3][ie])/(8.0*pow(pi, 2));
		Cp2_mixed += D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(Z_D[6][ie]*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je]) + 3.0*Z_D[6][je]*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie]))/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
		Cp3_mixed += D0mix*M_ch0[ae]*M_g*Z_D[6][ie]*pow(g_3, 2)*pow(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]), 2)*conj(Z_D[6][je])/(48.0*pow(pi, 2)) - D0mix*M_ch0[ae]*M_g*Z_D[6][je]*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*(conj(Z_N[1][ae])*sw/3.0 - conj(Z_N[2][ae])*cw)*conj(Z_D[gen][ie])/(2.0*cw*sw) + conj(Yd[gen])*conj(Z_D[gen+3][ie])*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])/(48.0*pow(pi, 2)) + D0mix*M_ch0[ae]*M_g*pow(g_3, 2)*(-sqrt(2)*Q_e*Z_D[6][ie]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][ie]*conj(Z_N[3][ae]))*(-sqrt(2)*Q_e*Z_D[6][je]*conj(Z_N[1][ae])/(3.0*cw) + Yd[3]*Z_D[3][je]*conj(Z_N[3][ae]))*conj(Z_D[gen][ie])*conj(Z_D[gen][je])/(48.0*pow(pi, 2));
    }
	
	CM_mixed[1]=C1_mixed;
	CM_mixed[2]=C2_mixed;
	CM_mixed[3]=C3_mixed;
	CM_mixed[4]=C4_mixed;
	CM_mixed[5]=C5_mixed;
	CMp_mixed[1]=Cp1_mixed;
	CMp_mixed[2]=Cp2_mixed;
	CMp_mixed[3]=Cp3_mixed;

    //higgs pinguin ->

    int ie,je;
	for(ie=1;ie<=5;ie++) CM_higgspenguin[ie];
	for(ie=1;ie<=3;ie++) CMp_higgspenguin[ie]=0.;

	double M_D[7],M_U[7];
    double complex Z_D[7][7],Z_U[7][7];

	M_D[1]=param->mass_dnl;
	M_D[2]=param->mass_stl;
	M_D[3]=param->mass_b1;
	M_D[4]=param->mass_dnr;
	M_D[5]=param->mass_str;
	M_D[6]=param->mass_b2;
	M_U[1]=param->mass_upl;
	M_U[2]=param->mass_chl;
	M_U[3]=param->mass_t1;
	M_U[4]=param->mass_upr;
	M_U[5]=param->mass_chr;
	M_U[6]=param->mass_t2;

	//In case sD_mix is not provided by the SLHA (set to 0):
	//getZD(Z_D,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_D[ie][je]= param->sD_mix[je][ie];
	//In case sU_mix is not provided by the SLHA (set to 0):
	//getZU(Z_U,param);
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) Z_U[ie][je]= conj(param->sU_mix[je][ie]);

	double complex delta_d[7][7],delta_d_LL[4][4],delta_d_LR[4][4],delta_d_RL[4][4],delta_d_RR[4][4];
	double complex delta_u[7][7],delta_u_LL[4][4],delta_u_LR[4][4],delta_u_RL[4][4],delta_u_RR[4][4];
	
	double m_av = (M_U[1]+M_U[2]+M_U[3]+M_U[4]+M_U[5]+M_U[6]+M_D[1]+M_D[2]+M_D[3]+M_D[4]+M_D[5]+M_D[6])/12.;
	
	getDelta(delta_d,Z_D,M_D,m_av,delta_d_LL,delta_d_LR,delta_d_RL,delta_d_RR);
	getDelta(delta_u,Z_U,M_U,m_av,delta_u_LL,delta_u_LR,delta_u_RL,delta_u_RR);
	
	double g_3=sqrt(4.*pi*alphas_running(mu_t,param->mass_top_pole,param->mass_b,param)); /* NM: compute from alphas instead of using g3 from SLHA file */
	double g_2=param->g2_Q;
	double M_g=param->mass_gluino;
	double tbeta=param->tan_beta;
	double M_W=param->mass_W;
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double m_t=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double M_A=param->mass_A0;
	double A_t=param->A_t;
	double mu = param->mu_Q;
	double complex M_2 = param->M2_Q;
	double x_mu = pow(abs(mu),2)/pow(m_av,2);
	double complex x_2 = pow(cabs(M_2),2)/pow(m_av,2);
	double x_g = pow(M_g,2)/pow(m_av,2);
	double alpha_s = pow(g_3,2.)/(4.*pi) ;
	double alpha_2 = pow(g_2,2.)/(4.*pi);
	double eps=2*alpha_s*mu*M_g*f(x_g)/(3.*pi*pow(m_av,2));
	
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];
	
	double complex V_tb = V_CKM[3][3]; 
	double complex V_tq = V_CKM[3][gen];
	double complex a1 = alpha_s*alpha_2*pow(m_b,2)*pow(tbeta,4)*pow(abs(mu),2)/(8*pi*pow(M_W,2)*pow(M_A,2)*pow(m_av,4)*pow((1+eps*tbeta),4));
	double complex a2 = -alpha_s*pow(M_g,2)*delta_d_LL[3][gen]*delta_d_RR[3][gen]*pow(h1(x_g),2);
	double complex a3 = alpha_2*pow(m_t,2)*A_t*M_g*h1(x_g)*h3(x_mu)*delta_d_RR[3][gen]*V_tb*conj(V_tq)/(pow(M_W,2));
	double complex a4 = alpha_2*M_2*M_g*delta_u_LL[3][gen]*delta_d_RR[3][gen]*h1(x_g)*h4(x_2,x_g);

	double complex C4_higgspenguin = a1*(a2+a3+a4);
	
	CM_higgspenguin[1]=0.;
	CM_higgspenguin[2]=0.;
	CM_higgspenguin[3]=0.;
	CM_higgspenguin[4]=C4_higgspenguin;
	CM_higgspenguin[5]=0.;
	CMp_higgspenguin[1]=0.;
	CMp_higgspenguin[2]=0.;
	CMp_higgspenguin[3]=0.;

}