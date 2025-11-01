	#include "include.h"

/*-------------------------------------------------------------------*/

double S0(double x) /* NM: needed for SM contribution */
{
	if(fabs(1.-x)<1.e-5) return S0(0.9999);

	return (4.*x-11.*x*x+x*x*x)/4./pow(1.-x,2.) - 3.*x*x*x*log(x)/2./pow(1.-x,3.);
}

/*-------------------------------------------------------------------*/

double D0(double w, double x, double y, double z)
{
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D0(0.9996,0.9998,1.0002,1.0004);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return D0(0.9996,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D0(w,0.9996,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D0(0.9996,x,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D0(0.9996,0.9998,y,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)) return D0(0.9998,1.0002,y,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)) return D0(0.9998,x,1.0002,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D0(0.9998,x,y,1.0002);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return D0(w,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D0(w,0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D0(w,x,0.9998,1.0002);
	if(fabs(1.-w)<1.e-5) return D0(0.9999,x,y,z);
	if(fabs(1.-x)<1.e-5) return D0(w,0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return D0(w,x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return D0(w,x,y,0.9999);
	if(fabs(1.-w/x)<1.e-5) return D0(x*0.9998,x,y,z);
	if(fabs(1.-w/y)<1.e-5) return D0(y*0.9998,x,y,z);
	if(fabs(1.-w/z)<1.e-5) return D0(z*0.9998,x,y,z);
	if(fabs(1.-x/y)<1.e-5) return D0(w,y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return D0(w,x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return D0(w,x,y,x*0.9998);

	return (w*log(w)/((z-w)*(y-w)*(x-w)))+(x*log(x)/((z-x)*(y-x)*(w-x)))+(y*log(y)/((z-y)*(w-y)*(x-y)))+(z*log(z)/((w-z)*(y-z)*(x-z)));
}

double D2p(double w, double x, double y, double z)
{
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D2p(0.9996,0.9998,1.0002,1.0004);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return D2p(0.9996,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D2p(w,0.9996,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D2p(0.9996,x,0.9998,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D2p(0.9996,0.9998,y,1.0002);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-x)<1.e-5)) return D2p(0.9998,1.0002,y,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-y)<1.e-5)) return D2p(0.9998,x,1.0002,z);
	if((fabs(1.-w)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D2p(0.9998,x,y,1.0002);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return D2p(w,0.9998,1.0002,z);
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D2p(w,0.9998,y,1.0002);
	if((fabs(1.-y)<1.e-5)&&(fabs(1.-z)<1.e-5)) return D2p(w,x,0.9998,1.0002);
	if(fabs(1.-w)<1.e-5) return D2p(0.9999,x,y,z);
	if(fabs(1.-x)<1.e-5) return D2p(w,0.9999,y,z);
	if(fabs(1.-y)<1.e-5) return D2p(w,x,0.9999,z);
	if(fabs(1.-z)<1.e-5) return D2p(w,x,y,0.9999);
	if(fabs(1.-w/x)<1.e-5) return D2p(x*0.9998,x,y,z);
	if(fabs(1.-w/y)<1.e-5) return D2p(y*0.9998,x,y,z);
	if(fabs(1.-w/z)<1.e-5) return D2p(z*0.9998,x,y,z);
	if(fabs(1.-x/y)<1.e-5) return D2p(w,y*0.9998,y,z);
	if(fabs(1.-y/z)<1.e-5) return D2p(w,x,z*0.9998,z);
	if(fabs(1.-x/z)<1.e-5) return D2p(w,x,y,x*0.9998);

	return 0.25*((w*w*log(w)/((z-w)*(y-w)*(x-w)))+(x*x*log(x)/((z-x)*(y-x)*(w-x)))+(y*y*log(y)/((z-y)*(w-y)*(x-y)))+(z*z*log(z)/((w-z)*(y-z)*(x-z))));
}

double h3(double x)
{
	if(fabs(x-1.)<1.e-5) return -1./4;
	return -0.5/(1.-x)-0.5*x*log(x)*pow(1.-x,-2.);
}

double h1(double x)
{
	if(fabs(x-1.)<1.e-5) return 4./9;
	return 4.*(1.+x)/(3.*pow(1.-x,2.)) + 8.*x*log(x)/(3.*pow((1.-x),3.));
}

double h4(double x, double y)
{
	if((fabs(1.-x)<1.e-5)&&(fabs(1.-y)<1.e-5)) return 1./6;
	if(fabs(1.-x)<1.e-5) return h4(0.9999,y);
	if(fabs(1.-y)<1.e-5) return h4(x,0.9999);
	if(fabs(1.-x/y)<1.e-5) return h4(y*0.9998,y);

	return -1./((1.-x)*(1.-y))+x*log(x)/((y-x)*pow((1.-x),2.))+y*log(y)/((x-y)*pow((1.-y),2.));
}

double f(double x)
{
	if(fabs(x-1.)<1.e-5) return 0.5;
	return 1./(1.-x)+x*log(x)/pow(1.-x,2.);
}

/*--------------------------------------------------------------------------*/

void getZD(double complex ZD[7][7], struct parameters* param)
{
	int ie,je;
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) ZD[ie][ie]=(ie==je);

	ZD[3][3] = param->sbot_mix[1][1];
	ZD[6][6] = ZD[3][3];
	ZD[3][6] = sin(acos(ZD[3][3]));
	ZD[6][3] = -ZD[3][6];
	
	return;
}

/*--------------------------------------------------------------------------*/

void getZU(double complex ZU[7][7], struct parameters* param)
{
	//double complex V_CKM[7][7]; /* NM: apparently no multiplication by CKM matrix */

	int ie,je,ke;
	//for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) V_CKM[ie][ie]=(ie==je);

    //for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];

	double complex ZU1[7][7];
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) ZU1[ie][ie]=(ie==je);

	ZU1[3][3] = param->stop_mix[1][1];
	ZU1[6][6] = ZU1[3][3];
	ZU1[3][6] = -sin(acos(ZU1[3][3]));
	ZU1[6][3] = -ZU1[3][6];

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		//ZU[ie][je] = 0.;
		//for(ke=1;ke<=6;ke++) ZU[ie][je]+=ZU1[ie][ke]*V_CKM[ke][je]; /* NM: apparently no multiplication by CKM matrix */
		for(ke=1;ke<=6;ke++) ZU[ie][je]=ZU1[ie][ke];
	}

	return;	
}

/*--------------------------------------------------------------------------*/

void getDelta(double complex delta[7][7],double complex Z[7][7],double M[7],double m_av,double complex delta_LL[4][4],double complex delta_LR[4][4],double complex delta_RL[4][4],double complex delta_RR[4][4])
{
	double complex mass_matrix[7][7];

	int ie,je,ke,le;
	
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		if (ie==je) mass_matrix[ie][je]=pow(M[ie],2.);
		else mass_matrix[ie][je]=0.;
	}
		
	double complex ZMZdag[7][7]; /* NM: ZMZdag = Z*M*Z^dagger */
	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++)
	{
		ZMZdag[ie][je]=0.;
		for(ke=1;ke<=6;ke++) for(le=1;le<=6;le++) ZMZdag[ie][je]+=Z[ie][ke]*mass_matrix[ke][le]*conj(Z[je][le]); /* NM: matrix product routine removed and replaced locally */
	}

	for(ie=1;ie<=6;ie++) for(je=1;je<=6;je++) delta[ie][je] = (ZMZdag[ie][je]-pow(m_av,2)*(ie==je))/(pow(m_av,2));
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++)
	{
		delta_LL[ie][je]=delta[ie][je];
		delta_LR[ie][je]=delta[ie][je+3];
		delta_RL[ie][je]=delta[ie+3][je];
		delta_RR[ie][je]=delta[ie+3][je+3];
	}
	
	return;
}

/*--------------------------------------------------------------------*/

void CM_calculator_sm(int gen, double complex CM_sm[], double complex CMp_sm[], double mu_t, struct parameters* param)

{
	int ie,je;
	for(ie=1;ie<=5;ie++) CM_sm[ie]=0.;
	for(ie=1;ie<=3;ie++) CMp_sm[ie]=0.;

	double complex V_tb = param->Vtb;
	double complex V_tq;

    V_tb = 0.999118;
	double G_F = param->Gfermi;
	double M_W = param->mass_W;
	double M_t = running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param);
	double m_B;
	
	double xt=pow(M_t/M_W,2.);
    xt = 4.663;
	if(gen==1) /* NM: added generation dependence */
	{
		V_tq = param->Vtd;
		m_B = param->m_Bd;
        V_tq = 0.00787861;
	}
	else
	{
		V_tq = param->Vts;
		m_B = param->m_Bs;
        V_tq = -0.0411012;
	}

	CM_sm[1]=pow(G_F*M_W,2.)*cpow(conj(V_tq)*V_tb,2.)*S0(xt)/(8.*pow(pi,2.)*m_B); /* NM: hardcoded S0 replaced by function */ 
	CM_sm[2]=0.;
	CM_sm[3]=0.;
	CM_sm[4]=0.;
	CM_sm[5]=0.;

	CMp_sm[1]=0.;
	CMp_sm[2]=0.;
	CMp_sm[3]=0.;

    printf("C_mix_bd_1 : %.4e\n", CM_sm[1]);
	return;
}

/* NM: Conventions for squark mixing matrices */
/* SLHA2: 0801.0045, see formulae (28) and (29) page 16 */
/* Buras: hep-ph/0703200, see formula (14) page 9 */
/* Apparently ZU_ij = conj(SLHAU_ji) and ZD_ij = SLHAD_ji */
/* Same matrices for neutralinos, but masses defined positive */
/* For charginos, Z+=V^dagger and Z-=Udagger */

void CM_calculator_gluino(int gen, double complex CM_gluino[], double complex CMp_gluino[], double mu_t, struct parameters* param)
{
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
	
	return;
}

void CM_calculator_chargedhiggs(int gen, double complex CM_chargedhiggs[],double complex CMp_chargedhiggs[], double mu_t, struct parameters* param)
{
	int ie,je;
	for(ie=1;ie<=5;ie++) CM_chargedhiggs[ie];
	for(ie=1;ie<=3;ie++) CMp_chargedhiggs[ie]=0.;

	double M_H=param->mass_H;
	double M_W=param->mass_W;
	double M_H_pow_2 = pow(M_H,2.);
	double M_W_pow_2 = pow(M_W,2.);
	
	double m_q;
	if(gen==1) m_q=param->mass_d;
	else m_q=param->mass_s;
	
	double m_b=running_mass(param->mass_b,param->mass_b,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	double g_2=param->g2_Q;	
	double tbeta=param->tan_beta;
	double m_u[4],m_u_pow_2[4];
	double complex V_CKM[4][4];
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++) V_CKM[ie][je]=param->CKM[ie][je]+I*param->IMCKM[ie][je];

	m_u[1]=param->mass_u;
	m_u[2]=running_mass(param->mass_c,param->mass_c,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */
	m_u[3]=running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param); /* NM: running mass */

	for(ie=1;ie<=3;ie++) m_u_pow_2[ie]=pow(m_u[ie],2.);
  
	double complex C1_chargedhiggs=0.;
	double complex C2_chargedhiggs=0.;
	double complex C3_chargedhiggs=0.;
	double complex C4_chargedhiggs=0.;
	double complex C5_chargedhiggs=0.;
	
	double complex Cp1_chargedhiggs=0.;
	double complex Cp2_chargedhiggs=0.;
	double complex Cp3_chargedhiggs=0.;
	
	double D0h,D0h_c,D2h,D2h_c;
	double complex CKM_product;
	
	for(ie=1;ie<=3;ie++) for(je=1;je<=3;je++)
	{
		D0h = D0(m_u_pow_2[ie],m_u_pow_2[je],M_H_pow_2,M_W_pow_2);
		D0h_c = D0(m_u_pow_2[ie],m_u_pow_2[je],M_H_pow_2,M_H_pow_2);
		D2h = D2p(m_u_pow_2[ie],m_u_pow_2[je],M_H_pow_2,M_W_pow_2);
		D2h_c = D2p(m_u_pow_2[ie],m_u_pow_2[je],M_H_pow_2,M_H_pow_2);
		
		CKM_product = V_CKM[ie][3]*V_CKM[je][3]*conj(V_CKM[ie][gen])*conj(V_CKM[je][gen]); 	/* NM: added generation dependence, V_CKM[ie][2] -> V_CKM[ie][gen] */

		C1_chargedhiggs += pow(g_2,4.)*CKM_product*m_u_pow_2[ie]*m_u_pow_2[je]*(2.*pow(M_W,2.)*D0h*pow(tbeta,-2.) - D2h_c*pow(tbeta,-4.) - 2*D2h*pow(tbeta,-2.))/(128.*pow(pi,2.)*pow(M_W,4.));
		C2_chargedhiggs += -pow(g_2,4.)*pow(m_q,2.)*CKM_product*pow(m_u[ie],2)*pow(m_u[je],2)*(D0h_c - 2*D0h)/(128.*pow(pi,2.)*pow(M_W,4.));
		C3_chargedhiggs += 0.;
		C4_chargedhiggs += pow(g_2,4.)*CKM_product*(m_b*m_q*D2h*pow(tbeta,2.)/pow(M_W,2.) - m_b*m_q*pow(m_u[ie],2)*pow(m_u[je],2)*(D0h_c + D0h*(pow(tbeta,2.) + pow(tbeta,-2.)))/(4*pow(M_W,4.)))/(16.*pow(pi,2.));
		C5_chargedhiggs += pow(g_2,4.)*m_b*m_q*CKM_product*pow(m_u[ie],2.)*(D2h_c-2.*D2h)/(32.*pow(pi,2.)*pow(M_W,4.));
		Cp1_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*pow(m_q,2.)*CKM_product*(D2h_c*pow(tbeta,4.)+2.*D2h*pow(tbeta,2.))/(128.*pow(pi,2.)*pow(M_W,4.));
		Cp2_chargedhiggs += -pow(g_2,4.)*pow(m_b,2.)*CKM_product*pow(m_u[ie],2)*pow(m_u[je],2)*(D0h_c-2.*D0h)/(128.*pow(pi,2.)*pow(M_W,4.));
		Cp3_chargedhiggs += 0.;
	}
	CM_chargedhiggs[1]=C1_chargedhiggs;
	CM_chargedhiggs[2]=C2_chargedhiggs;
	CM_chargedhiggs[3]=C3_chargedhiggs;
	CM_chargedhiggs[4]=C4_chargedhiggs;
	CM_chargedhiggs[5]=C5_chargedhiggs;
	
	CMp_chargedhiggs[1]=Cp1_chargedhiggs;
	CMp_chargedhiggs[2]=Cp2_chargedhiggs;
	CMp_chargedhiggs[3]=Cp3_chargedhiggs;

	return;
}

void CM_calculator_chargino(int gen, double complex CM_chargino[], double complex CMp_chargino[], double mu_t, struct parameters* param)
{
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
	
	return;
}

void CM_calculator_neutralino(int gen, double complex CM_neutralino[], double complex CMp_neutralino[], double mu_t, struct parameters* param)
{
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
	
	return;
}

void CM_calculator_mixed(int gen, double complex CM_mixed[], double complex CMp_mixed[], double mu_t, struct parameters* param)
{
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
	
	return;
}

void CM_calculator_higgspenguin(int gen, double complex CM_higgspenguin[], double complex CMp_higgspenguin[], double mu_t, struct parameters* param)
/* NM: where is it taken from? */
{
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
	return;
}

void CM_running(double complex CM[],double complex CMp[], double mu_t, double complex CMb[], double complex CMpb[], double mu_b, struct parameters* param)
{
	int ie,je;
	for(ie=1;ie<=5;ie++) CMb[ie]=CMpb[ie]=0.;

	double alphas_mut=alphas_running(mu_t,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mu=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);	
	double eta=alphas_mut/alphas_mu;
	
	/* LO from 1704.06639 */
	double M[8][8]={{pow(eta,6./23.), 0., 0., 0., 0., 0., 0., 0.},
		           {0., 0.983117/pow(eta,0.631486) + 0.0168825*pow(eta,0.718442), -(0.257663/pow(eta,0.631486)) + 0.257663*pow(eta,0.718442), 0., 0., 0., 0., 0.},
		           {0., -(0.0644157/pow(eta,0.631486)) + 0.0644157*pow(eta,0.718442), 0.0168825/pow(eta,0.631486) + 0.983117*pow(eta,0.718442), 0., 0., 0., 0., 0.},
		           {0., 0., 0., 1./pow(eta,24./23.), 0.333333/pow(eta,24./23.) - 0.333333*pow(eta,3./23.), 0., 0., 0.},
		           {0., 0., 0., 0., pow(eta,3./23.), 0., 0., 0.},
		           {0., 0., 0., 0., 0., pow(eta,6./23.), 0., 0.},
		           {0., 0., 0., 0., 0., 0., 0.983117/pow(eta,0.631486) + 0.0168825*pow(eta,0.718442), -(0.257663/pow(eta,0.631486)) + 0.257663*pow(eta,0.718442)},
		           {0., 0., 0., 0., 0., 0., -(0.0644157/pow(eta,0.631486)) + 0.0644157*pow(eta,0.718442), 0.0168825/pow(eta,0.631486) + 0.983117*pow(eta,0.718442)}};



	double complex C[8]={CM[1],CM[2],CM[3],CM[4],CM[5],CMp[1],CMp[2],CMp[3]};

	double complex Cb[8];
	for(ie=0;ie<8;ie++) Cb[ie]=0.;
	for(ie=0;ie<8;ie++) for(je=0;je<8;je++) Cb[ie]+=M[ie][je]*C[je];

	for(ie=1;ie<=5;ie++) CMb[ie]=Cb[ie-1];
	for(ie=1;ie<=3;ie++) CMpb[ie]=Cb[ie+4];
	
	CMb[1]+=alphas_mu/4./pi*(1.6273*(1.-eta)*pow(eta,6./23.))*CM[1]; /* NLO from hep-ph/0102316 */ //SN: seems to be eq.3.1 but only for VLL and the rest remain at LO
	CMpb[1]+=alphas_mu/4./pi*(1.6273*(1.-eta)*pow(eta,6./23.))*CMp[1];

	return;
}

void CM_calculator(int gen, double complex CM[],double complex CMp[], double mu_t, struct parameters* param, int model, int gluino, int chargino, int charghiggs, int mixed, int neutralino, int penguin)
{
	int ie;

	double complex CMt[6];
	double complex CMpt[4];
	for(ie=1;ie<=5;ie++) CMt[ie]=0.;
	for(ie=1;ie<=3;ie++) CMpt[ie]=0.;
	
	/* WC formulae from hep-ph/0703200 */
	double complex CM_gluino[6];
	double complex CMp_gluino[4];
	
	double complex CM_chargedhiggs[6];
	double complex CMp_chargedhiggs[4];
	
	double complex CM_chargino[6];
	double complex CMp_chargino[4];

	double complex CM_neutralino[6];
	double complex CMp_neutralino[4];
	
	double complex CM_mixed[6];
	double complex CMp_mixed[4];
	
	double complex CM_higgspenguin[6];
	double complex CMp_higgspenguin[4];

	double complex CM_sm[6];
	double complex CMp_sm[4];

	double mB;
	if(gen==1) mB=param->m_Bd; else mB=param->m_Bs;

	CM_calculator_sm(gen,CM_sm,CMp_sm,mu_t,param);
	
	if(model==1) /* NM: MSSM */
	{
		CM_calculator_gluino(gen,CM_gluino,CMp_gluino,mu_t,param);
		CM_calculator_chargedhiggs(gen,CM_chargedhiggs,CMp_chargedhiggs,mu_t,param);
		CM_calculator_chargino(gen,CM_chargino,CMp_chargino,mu_t,param);
		CM_calculator_neutralino(gen,CM_neutralino,CMp_neutralino,mu_t,param);
		CM_calculator_mixed(gen,CM_mixed,CMp_mixed,mu_t,param);
		CM_calculator_higgspenguin(gen,CM_higgspenguin,CMp_higgspenguin,mu_t,param);

		CM[1]=CM_sm[1]+(gluino*CM_gluino[1]+mixed*CM_mixed[1]+neutralino*CM_neutralino[1]+chargino*CM_chargino[1]+charghiggs*CM_chargedhiggs[1]+penguin*CM_higgspenguin[1])/(2.*mB);
		CM[2]=CM_sm[2]+(gluino*CM_gluino[2]+mixed*CM_mixed[2]+neutralino*CM_neutralino[2]+chargino*CM_chargino[2]+charghiggs*CM_chargedhiggs[2]+penguin*CM_higgspenguin[2])/(2.*mB);
		CM[3]=CM_sm[3]+(gluino*CM_gluino[3]+mixed*CM_mixed[3]+neutralino*CM_neutralino[3]+chargino*CM_chargino[3]+charghiggs*CM_chargedhiggs[3]+penguin*CM_higgspenguin[3])/(2.*mB);
		CM[4]=CM_sm[4]+(gluino*CM_gluino[4]+mixed*CM_mixed[4]+neutralino*CM_neutralino[4]+chargino*CM_chargino[4]+charghiggs*CM_chargedhiggs[4]+penguin*CM_higgspenguin[4])/(2.*mB);
		CM[5]=CM_sm[5]+(gluino*CM_gluino[5]+mixed*CM_mixed[5]+neutralino*CM_neutralino[5]+chargino*CM_chargino[5]+charghiggs*CM_chargedhiggs[5]+penguin*CM_higgspenguin[5])/(2.*mB);

		CMp[1]=CMp_sm[1]+(gluino*CMp_gluino[1]+mixed*CMp_mixed[1]+neutralino*CMp_neutralino[1]+chargino*CMp_chargino[1]+charghiggs*CMp_chargedhiggs[1]+penguin*CMp_higgspenguin[1])/(2.*mB);
		CMp[2]=CMp_sm[2]+(gluino*CMp_gluino[2]+mixed*CMp_mixed[2]+neutralino*CMp_neutralino[2]+chargino*CMp_chargino[2]+charghiggs*CMp_chargedhiggs[2]+penguin*CMp_higgspenguin[2])/(2.*mB);
		CMp[3]=CMp_sm[3]+(gluino*CMp_gluino[3]+mixed*CMp_mixed[3]+neutralino*CMp_neutralino[3]+chargino*CMp_chargino[3]+charghiggs*CMp_chargedhiggs[3]+penguin*CMp_higgspenguin[3])/(2.*mB);		
	}
	else if(model==0) /* NM: SM only */
	{
		CM[1]=CM_sm[1];
		CM[2]=CM_sm[2];
		CM[3]=CM_sm[3];
		CM[4]=CM_sm[4];
		CM[5]=CM_sm[5];

		CMp[1]=CMp_sm[1];
		CMp[2]=CMp_sm[2];
		CMp[3]=CMp_sm[3];
	}
  	else if(model==-1) /* Wilson coefficients input */
	{
		// CM[1]=CM_sm[1] + (param->CM_Bmix[gen][1]); //THEO : no input
		// CM[2]=CM_sm[2]+ (param->CM_Bmix[gen][2]);
		// CM[3]=CM_sm[3]+ (param->CM_Bmix[gen][3]);
		// CM[4]=CM_sm[4]+ (param->CM_Bmix[gen][4]);
		// CM[5]=CM_sm[5]+ (param->CM_Bmix[gen][5]);

		// CMp[1]=CMp_sm[1]+ (param->CMp_Bmix[gen][1]);
		// CMp[2]=CMp_sm[2]+ (param->CMp_Bmix[gen][2]);
		// CMp[3]=CMp_sm[3]+ (param->CMp_Bmix[gen][3]);

		// printf("Using CM1 is: %.14e %.14e\n",CM_sm[1] , (param->CM_Bmix[gen][1]));
	}
	return;
}

/*********************************************************/
/*--------------------------------------------------------------------*/


void CM_calculatorB(int gen, double complex CM[],double complex CMp[], double mu_t, struct parameters* param, int model, int gluino, int chargino, int charghiggs, int mixed, int neutralino, int penguin)
{
	int ie;

	double complex CMt[6];
	double complex CMpt[4];
	for(ie=1;ie<=5;ie++) CMt[ie]=0.;
	for(ie=1;ie<=3;ie++) CMpt[ie]=0.;

	/* WC formulae from hep-ph/0703200 */
	double complex CM_gluino[6];
	double complex CMp_gluino[4];

	double complex CM_chargedhiggs[6];
	double complex CMp_chargedhiggs[4];

	double complex CM_chargino[6];
	double complex CMp_chargino[4];

	double complex CM_neutralino[6];
	double complex CMp_neutralino[4];

	double complex CM_mixed[6];
	double complex CMp_mixed[4];

	double complex CM_higgspenguin[6];
	double complex CMp_higgspenguin[4];

	double complex CM_sm[6];
	double complex CMp_sm[4];

	double mB;
	if(gen==1) mB=param->m_Bd; else mB=param->m_Bs;

	/*SN: putting SM WC to zero; M12 for SM is directly cacluated with M12_calculatorB_sm in mixing.c*/
	// CM_calculator_sm(gen,CM_sm,CMp_sm,mu_t,param);
	for(ie=0;ie<6;ie++) CM_sm[ie] = 0.;
	for(ie=0;ie<4;ie++) CMp_sm[ie] = 0.;

	if(model==1) /* NM: MSSM */
	{
		CM_calculator_gluino(gen,CM_gluino,CMp_gluino,mu_t,param);
		CM_calculator_chargedhiggs(gen,CM_chargedhiggs,CMp_chargedhiggs,mu_t,param);
		CM_calculator_chargino(gen,CM_chargino,CMp_chargino,mu_t,param);
		CM_calculator_neutralino(gen,CM_neutralino,CMp_neutralino,mu_t,param);
		CM_calculator_mixed(gen,CM_mixed,CMp_mixed,mu_t,param);
		CM_calculator_higgspenguin(gen,CM_higgspenguin,CMp_higgspenguin,mu_t,param);

		CM[1]=CM_sm[1]+(gluino*CM_gluino[1]+mixed*CM_mixed[1]+neutralino*CM_neutralino[1]+chargino*CM_chargino[1]+charghiggs*CM_chargedhiggs[1]+penguin*CM_higgspenguin[1])/(2.*mB);
		CM[2]=CM_sm[2]+(gluino*CM_gluino[2]+mixed*CM_mixed[2]+neutralino*CM_neutralino[2]+chargino*CM_chargino[2]+charghiggs*CM_chargedhiggs[2]+penguin*CM_higgspenguin[2])/(2.*mB);
		CM[3]=CM_sm[3]+(gluino*CM_gluino[3]+mixed*CM_mixed[3]+neutralino*CM_neutralino[3]+chargino*CM_chargino[3]+charghiggs*CM_chargedhiggs[3]+penguin*CM_higgspenguin[3])/(2.*mB);
		CM[4]=CM_sm[4]+(gluino*CM_gluino[4]+mixed*CM_mixed[4]+neutralino*CM_neutralino[4]+chargino*CM_chargino[4]+charghiggs*CM_chargedhiggs[4]+penguin*CM_higgspenguin[4])/(2.*mB);
		CM[5]=CM_sm[5]+(gluino*CM_gluino[5]+mixed*CM_mixed[5]+neutralino*CM_neutralino[5]+chargino*CM_chargino[5]+charghiggs*CM_chargedhiggs[5]+penguin*CM_higgspenguin[5])/(2.*mB);

		CMp[1]=CMp_sm[1]+(gluino*CMp_gluino[1]+mixed*CMp_mixed[1]+neutralino*CMp_neutralino[1]+chargino*CMp_chargino[1]+charghiggs*CMp_chargedhiggs[1]+penguin*CMp_higgspenguin[1])/(2.*mB);
		CMp[2]=CMp_sm[2]+(gluino*CMp_gluino[2]+mixed*CMp_mixed[2]+neutralino*CMp_neutralino[2]+chargino*CMp_chargino[2]+charghiggs*CMp_chargedhiggs[2]+penguin*CMp_higgspenguin[2])/(2.*mB);
		CMp[3]=CMp_sm[3]+(gluino*CMp_gluino[3]+mixed*CMp_mixed[3]+neutralino*CMp_neutralino[3]+chargino*CMp_chargino[3]+charghiggs*CMp_chargedhiggs[3]+penguin*CMp_higgspenguin[3])/(2.*mB);
	}
	else if(model==0)
	{
		CM[1]=CM_sm[1];
		CM[2]=CM_sm[2];
		CM[3]=CM_sm[3];
		CM[4]=CM_sm[4];
		CM[5]=CM_sm[5];

		CMp[1]=CMp_sm[1];
		CMp[2]=CMp_sm[2];
		CMp[3]=CMp_sm[3];
	}
	else if(model==-1) /* Wilson coefficients input */
	{
		// CM[1]=CM_sm[1] + (param->CM_Bmix[gen][1]);;
		// CM[2]=CM_sm[2]+ (param->CM_Bmix[gen][2]);; //THEO : no input
		// CM[3]=CM_sm[3]+ (param->CM_Bmix[gen][3]);;
		// CM[4]=CM_sm[4]+ (param->CM_Bmix[gen][4]);;
		// CM[5]=CM_sm[5]+ (param->CM_Bmix[gen][5]);;

		// CMp[1]=CMp_sm[1]+ (param->CMp_Bmix[gen][1]);;
		// CMp[2]=CMp_sm[2]+ (param->CMp_Bmix[gen][2]);;
		// CMp[3]=CMp_sm[3]+ (param->CMp_Bmix[gen][3]);;
	}
	return;
}



void CM_calculatorK_sm(double complex CM_sm[], double complex CMp_sm[], double mu_t, struct parameters* param)

{
	int ie,je;
	for(ie=1;ie<=5;ie++) CM_sm[ie]=0;
	for(ie=1;ie<=3;ie++) CMp_sm[ie]=0.;

	double complex V_td = param->Vtd;
	double complex V_ts = param->Vts;
	double G_F = param->Gfermi;
	double M_W = param->mass_W;
	double M_t = running_mass(param->mtmt,param->mtmt,mu_t,param->mass_top_pole,param->mass_b,param);
	double m_K = param->m_K0;

	double xt=pow(M_t/M_W,2.);

	/*SN: TO BE CHANGED */
	CM_sm[1]=0.;//pow(G_F*M_W,2.)*cpow(conj(V_td)*V_ts,2.)*S0(xt)/(8.*pow(pi,2.)*m_K);
	CM_sm[2]=0.;
	CM_sm[3]=0.;
	CM_sm[4]=0.;
	CM_sm[5]=0.;

	CMp_sm[1]=0.;
	CMp_sm[2]=0.;
	CMp_sm[3]=0.;

	return;
}




void CM_calculatorK(double complex CM[],double complex CMp[], double mu_t, struct parameters* param, int model)
{
/* The Wilson coefficients here are given in SUSY basis */

	double complex CM_sm[6];
	double complex CMp_sm[4];

	CM_calculatorK_sm(CM_sm,CMp_sm,mu_t,param);

	if(model==1) /* SN: general */
	{
		CM[1]= CM_sm[1] + 0.;
		CM[2]= 0.;
		CM[3]= 0.;
		CM[4]= 0.;
		CM[5]= 0.;

		CMp[1]= 0.;
		CMp[2]= 0.;
		CMp[3]= 0.;
	}
	else if(model==0) /* SM only */
	{
		CM[1]=CM_sm[1];
		CM[2]=CM_sm[2];
		CM[3]=CM_sm[3];
		CM[4]=CM_sm[4];
		CM[5]=CM_sm[5];

		CMp[1]=CMp_sm[1];
		CMp[2]=CMp_sm[2];
		CMp[3]=CMp_sm[3];
	}
	else if(model==-1) /* Wilson coefficients input */
	{
		// CM[1]=CM_sm[1] + (param->CM_Kmix[1]);  //THEO : no input
		// CM[2]=CM_sm[2]+ (param->CM_Kmix[2]);
		// CM[3]=CM_sm[3]+ (param->CM_Kmix[3]);
		// CM[4]=CM_sm[4]+ (param->CM_Kmix[4]);
		// CM[5]=CM_sm[5]+ (param->CM_Kmix[5]);

		// CMp[1]=CMp_sm[1]+ (param->CMp_Kmix[1]);
		// CMp[2]=CMp_sm[2]+ (param->CMp_Kmix[2]);
		// CMp[3]=CMp_sm[3]+ (param->CMp_Kmix[3]);

		
	}
	return;
}

void CM_calculatorD(double complex CM[],double complex CMp[], double mu_t, struct parameters* param, int model)
{
/* The Wilson coefficients here are given in SUSY basis */

	double complex CM_sm[6];
	double complex CMp_sm[4];

	int i;
	for(i=0;i<=5;i++) CM_sm[i]=0.;
	for(i=0;i<=3;i++) CMp_sm[i]=0.;

	if(model==1) /* SN: general */
	{
		CM[1]= 0.;
		CM[2]= 0.;
		CM[3]= 0.;
		CM[4]= 0.;
		CM[5]= 0.;

		CMp[1]= 0.;
		CMp[2]= 0.;
		CMp[3]= 0.;
	}
	else if(model==0) /* SM only */
	{
		CM[1]=CM_sm[1];
		CM[2]=CM_sm[2];
		CM[3]=CM_sm[3];
		CM[4]=CM_sm[4];
		CM[5]=CM_sm[5];

		CMp[1]=CMp_sm[1];
		CMp[2]=CMp_sm[2];
		CMp[3]=CMp_sm[3];
	}
	else if(model==-1)  /* Wilson coefficients input */
	{
		// CM[1]=CM_sm[1] + (param->CM_Dmix[1]);  //THEO : no input
		// CM[2]=CM_sm[2] + (param->CM_Dmix[2]);
		// CM[3]=CM_sm[3] + (param->CM_Dmix[3]);
		// CM[4]=CM_sm[4] + (param->CM_Dmix[4]);
		// CM[5]=CM_sm[5] + (param->CM_Dmix[5]);

		// CMp[1]=CMp_sm[1] + (param->CMp_Dmix[1]);
		// CMp[2]=CMp_sm[2] + (param->CMp_Dmix[2]);
		// CMp[3]=CMp_sm[3] + (param->CMp_Dmix[3]);

		// printf("Using CM1 D mixing is: %.14e %.14e\n",cabs(CM_sm[1]) , cabs(param->CM_Dmix[1]));
	}
	return;
}

void CM_converter_SUSYtoBMU(double complex CM[],double complex CMp[],double complex CM_BMU[])
{
	// VLL->Q[0], SLR1->Q[1], SLR2->Q[2], SLL1->Q[3], SLL2->Q[4], VRR->Q[5], SRR1->Q[6], SRR2->Q[7]
	// double complex CM_BMU[9];
/*VLL*/		CM_BMU[0] = CM[1];
/*LR1*/		CM_BMU[1] = -1./2.*CM[5];
/*LR2*/		CM_BMU[2] = CM[4];
/*SLL1*/	CM_BMU[3] = CM[2]-1./2.*CM[3];//or perhaps = CM[2]+5./2.*CM[3]
/*SLL2*/	CM_BMU[4] = 1./8.*CM[3];//or perhaps = -3./2.*CM[3]
/*VRR*/		CM_BMU[5] = CMp[1];
/*SRR1*/	CM_BMU[6] = CMp[2]-1./2.*CMp[3];//or perhaps = CMp[2]+5./2.*CMp[3]
/*SRR2*/	CM_BMU[7] = 1./8.*CMp[3];//or perhaps = -3./2.*CMp[3]

	return;
}


void CM_running_tb_BMU(double complex CMt_BMU[], double mu_t, double complex CMb_BMU[], double mu_b, struct parameters* param)
{
	double alphas_mut=alphas_running(mu_t,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	double eta5=alphas_mut/alphas_mub;
	printf("as_ew=%.5f, as_b=%.5f\n", alphas_mut, alphas_mub);

/************************************************************/
/* LO running from mu_t to mu_b in the BMU basis 0102316 */
	double M_tb_LO_BMU[8][8]={{pow(eta5,6./23.), 0., 0., 0., 0., 0., 0., 0.},
		                    {0., pow(eta5,3./23.), 0., 0., 0., 0., 0., 0.},
		                    {0., 2./3.*(pow(eta5,3./23.)-pow(eta5,-24./23.)), pow(eta5,-24./23.), 0., 0., 0., 0., 0.},
		                    {0., 0., 0., 1.0153*pow(eta5,-0.6315)-0.0153*pow(eta5,0.7184), 1.9325*(pow(eta5,-0.6315)-pow(eta5,0.7184)), 0., 0., 0.},
		                    {0., 0., 0., 0.0081*(pow(eta5,0.7184)-pow(eta5,-0.6315)), 1.0153*pow(eta5,0.7184)-0.0153*pow(eta5,-0.6315), 0., 0., 0.},
		                    {0., 0., 0., 0., 0., pow(eta5,6./23.), 0., 0.},
		                    {0., 0., 0., 0., 0., 0., 1.0153*pow(eta5,-0.6315)-0.0153*pow(eta5,0.7184), 1.9325*(pow(eta5,-0.6315)-pow(eta5,0.7184))},
		                    {0., 0., 0., 0., 0., 0., 0.0081*(pow(eta5,0.7184)-pow(eta5,-0.6315)), 1.0153*pow(eta5,0.7184)-0.0153*pow(eta5,-0.6315)}};
/***********************************/
/* NLO running from mu_t to mu_b in the BMU basis 0102316 */
	double M_tb_NLO_BMU[8][8]={{1.6273*(1.-eta5)*pow(eta5,6./23.), 0., 0., 0., 0., 0., 0., 0.},
		                     {0., 0.9250*pow(eta5,-24./23.)+pow(eta5,3./23.)*(-2.0994+1.1744*eta5), 1.3875*(pow(eta5,26./23.)-pow(eta5,-24./23.)), 0., 0., 0., 0., 0.},
		                     {0., (-11.7329+0.7829*eta5)*pow(eta5,3./23.)+pow(eta5,-24./23.)*(-5.3048+16.2548*eta5), (7.9572-8.8822*eta5)*pow(eta5,-24./23.)+0.9250*pow(eta5,26./23.), 0., 0., 0., 0., 0.},
		                     {0., 0., 0., (4.8177-5.2272*eta5)*pow(eta5,-0.6315)+(0.3371+0.0724*eta5)*pow(eta5,0.7184), (9.1696-38.8778*eta5)*pow(eta5,-0.6315)+(42.5021-12.7939*eta5)*pow(eta5,0.7184), 0., 0., 0.},
		                     {0., 0., 0., (0.0531+0.0415*eta5)*pow(eta5,-0.6315)-(0.0566 + 0.0380*eta5)*pow(eta5,0.7184), (0.1011+0.3083*eta5)*pow(eta5,-0.6315)+(-7.1314 + 6.7219*eta5)*pow(eta5,0.7184), 0., 0., 0.},
		                     {0., 0., 0., 0., 0., 1.6273*(1.-eta5)*pow(eta5,6./23.), 0., 0.},
		                     {0., 0., 0., 0., 0., 0., (4.8177-5.2272*eta5)*pow(eta5,-0.6315)+(0.3371+0.0724*eta5)*pow(eta5,0.7184), (9.1696-38.8778*eta5)*pow(eta5,-0.6315)+(42.5021-12.7939*eta5)*pow(eta5,0.7184)},
		                     {0., 0., 0., 0., 0., 0., (0.0531+0.0415*eta5)*pow(eta5,-0.6315)-(0.0566 + 0.0380*eta5)*pow(eta5,0.7184), (0.1011+0.3083*eta5)*pow(eta5,-0.6315)+(-7.1314 + 6.7219*eta5)*pow(eta5,0.7184)}};
/************************************************************/

	int ie,je;
	for(ie=0;ie<8;ie++) CMb_BMU[ie]=0.;
	for(ie=0;ie<8;ie++) for(je=0;je<8;je++) CMb_BMU[ie]+=(M_tb_LO_BMU[ie][je]+alphas_mub/4./pi*M_tb_NLO_BMU[ie][je])*CMt_BMU[je];

	return;
}

void CM_running_tL_BMU(double complex CMt_BMU[], double mu_t, double complex CML_BMU[], double mu_b, double mu_L, struct parameters* param)
{
	double alphas_mut=alphas_running(mu_t,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mub=alphas_running(mu_b,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_muL=alphas_running(mu_L,param->mass_top_pole,param->mass_b_pole,param);
	double eta5=alphas_mut/alphas_mub;
	double eta4=alphas_mub/alphas_muL;
	printf("as_ew=%.5f, as_b=%.5f, as_L=%.5f\n", alphas_mut, alphas_mub, alphas_muL);

/******************************************************************/
/* LO running from mu_t to mu_L ~ 2-3 GeV in the BMU basis 0102316 */

/*
*  i		1		2		3		4			5
*  a_i		6/23	3/23	-24/23	-0.6315		0.7184
*  b_i		6/25	3/25	-24/25	-0.5810		0.6610	
*/

double e41 = pow(eta4,6./25.);
double e42 = pow(eta4,3./25.);
double e43 = pow(eta4,-24./25.);
double e44 = pow(eta4,-0.5810);
double e45 = pow(eta4,0.6610);
double e51 = pow(eta5,6./23.);
double e52 = pow(eta5,3./23.);
double e53 = pow(eta5,-24./23.);
double e54 = pow(eta5,-0.6315);
double e55 = pow(eta5,0.7184);

	double M_tL_LO_BMU[8][8]={
		{e41*e51, 0.		 			 , 0.     , 0.							 , 0.							, 0.	 , 0.							, 0.						   },
		{0.		, e42*e52				 , 0.     , 0.							 , 0.							, 0.	 , 0.							, 0.						   },
		{0.		, 2./3.*(e42*e52-e43*e53), e43*e53, 0.							 , 0.							, 0.	 , 0.							, 0.						   },
		{0.		, 0.					 , 0.     , 1.0153*e44*e54-0.0153*e45*e55, 1.9325*(e44*e54-e45*e55)		, 0.	 , 0.							, 0.						   },
		{0.		, 0.					 , 0.	  , 0.0081*(e45*e55-e44*e54)	 , 1.0153*e45*e55-0.0153*e44*e54, 0.	 , 0.							, 0.						   },
		{0.		, 0.					 , 0.	  , 0.							 , 0.							, e41*e51, 0.							, 0.						   },
		{0.		, 0.					 , 0.	  , 0.							 , 0.							, 0.	 , 1.0153*e44*e54-0.0153*e45*e55, 1.9325*(e44*e54-e45*e55)	   },
		{0.		, 0.					 , 0.	  , 0.							 , 0.							, 0.	 , 0.0081*(e45*e55-e44*e54)		, 1.0153*e45*e55-0.0153*e44*e54}};
/***********************************/
/* NLO running from mu_t to mu_L ~ 2-3 GeV in the BMU basis 0102316 */
	double M_tL_NLO_BMU[8][8]={{pow(eta4,6./25.)*pow(eta5,6./23.)*(1.7917-0.1644*eta4-1.6273*eta4*eta5), 0., 0., 0., 0., 0., 0., 0.},
		                    {0., 0.9279*pow(eta4,-24./25.)*pow(eta5,-24./23.)-0.0029*pow(eta4,28./25.)*pow(eta5,-24./23.)+pow(eta4,3./25.)*pow(eta5,3./23.)*(-2.0241-0.0753*eta4+1.1744*eta4*eta5), -1.3918*pow(eta4,-24./25.)*pow(eta5,-24./23.)+0.0043*pow(eta4,28./25.)*pow(eta5,-24./23.)+1.3875*pow(eta4,28./25.)*pow(eta5,26./23.), 0., 0., 0., 0., 0.},
		                    {0., -0.0019*pow(eta4,28./25.)*pow(eta5,-24./23.) + 5.0000*pow(eta4,1./25.)*pow(eta5,3./23.)+pow(eta4,3./25.)*pow(eta5,3./23.)*(-16.6828-0.0502*eta4+0.7829*eta4*eta5)+pow(eta4,-24./25.)*pow(eta5,-24./23.)*(-4.4701-0.8327*eta4+16.2548*eta4*eta5 ), 0.0029*pow(eta4,28./25.)*pow(eta5,-24./23.)+0.9250*pow(eta4,28./25.)*pow(eta5,26./23.)+pow(eta4,-24./25.)*pow(eta5,-24./23.)*(6.7052+1.2491*eta4-8.8822*eta4*eta5), 0., 0., 0., 0., 0.},
		                    {0., 0., 0., 0.0020*pow(eta4,1.6610)*pow(eta5,-0.6315)-0.0334*pow(eta4,0.4190)*pow(eta5,0.7184)+pow(eta4,-0.5810)*pow(eta5,-0.6315)*(4.2458+0.5700*eta4-5.2272*eta4*eta5)+pow(eta4,0.6610)*pow(eta5,0.7184)*(0.3640+0.0064*eta4+0.0724*eta4*eta5), 0.0038*pow(eta4,1.6610)*pow(eta5,-0.6315)-4.2075*pow(eta4, 0.4190)*pow(eta5,0.7184)+pow(eta4, -0.5810)*pow(eta5,-0.6315)*(8.0810 + 1.0848*eta4-38.8778*eta4*eta5)+pow(eta4,0.6610)*pow(eta5,0.7184)*(45.9008+0.8087*eta4-12.7939*eta4*eta5), 0., 0., 0.},
		                    {0., 0., 0., -0.0011*pow(eta4,1.6610)*pow(eta5,-0.6315)+0.0003*pow(eta4,0.4190)*pow(eta5,0.7184)+pow(eta4,0.6610)*pow(eta5,0.7184)*(-0.0534-0.0034*eta4-0.0380*eta4*eta5)+pow(eta4,-0.5810)*pow(eta5,-0.6315)*(0.0587-0.0045*eta4+0.0415*eta4*eta5 ), -0.0020*pow(eta4,1.6610)*pow(eta5,-0.6315)+0.0334*pow(eta4,0.4190)*pow(eta5,0.7184)+pow(eta4,-0.5810)*pow(eta5,-0.6315)*(0.1117-0.0086*eta4+0.3083*eta4*eta5)+pow(eta4,0.6610)*pow(eta5,0.7184)*(-6.7398-0.4249*eta4+6.7219*eta4*eta5), 0., 0., 0.},
		                    {0., 0., 0., 0., 0., pow(eta4,6./25.)*pow(eta5,6./23.)*(1.7917-0.1644*eta4-1.6273*eta4*eta5), 0., 0.},
		                    {0., 0., 0., 0., 0., 0., 0.0020*pow(eta4,1.6610)*pow(eta5,-0.6315)-0.0334*pow(eta4,0.4190)*pow(eta5,0.7184)+pow(eta4,-0.5810)*pow(eta5,-0.6315)*(4.2458+0.5700*eta4-5.2272*eta4*eta5)+pow(eta4,0.6610)*pow(eta5,0.7184)*(0.3640+0.0064*eta4+0.0724*eta4*eta5), 0.0038*pow(eta4,1.6610)*pow(eta5,-0.6315)-4.2075*pow(eta4, 0.4190)*pow(eta5,0.7184)+pow(eta4, -0.5810)*pow(eta5,-0.6315)*(8.0810 + 1.0848*eta4-38.8778*eta4*eta5)+pow(eta4,0.6610)*pow(eta5,0.7184)*(45.9008+0.8087*eta4-12.7939*eta4*eta5)},
		                    {0., 0., 0., 0., 0., 0., -0.0011*pow(eta4,1.6610)*pow(eta5,-0.6315)+0.0003*pow(eta4,0.4190)*pow(eta5,0.7184)+pow(eta4,0.6610)*pow(eta5,0.7184)*(-0.0534-0.0034*eta4-0.0380*eta4*eta5)+pow(eta4,-0.5810)*pow(eta5,-0.6315)*(0.0587-0.0045*eta4+0.0415*eta4*eta5 ), -0.0020*pow(eta4,1.6610)*pow(eta5,-0.6315)+0.0334*pow(eta4,0.4190)*pow(eta5,0.7184)+pow(eta4,-0.5810)*pow(eta5,-0.6315)*(0.1117-0.0086*eta4+0.3083*eta4*eta5)+pow(eta4,0.6610)*pow(eta5,0.7184)*(-6.7398-0.4249*eta4+6.7219*eta4*eta5)}};
/******************************************************************/

	int ie,je;
	for(ie=0;ie<8;ie++) CML_BMU[ie]=0.;

	for(ie=0;ie<8;ie++) for(je=0;je<8;je++) CML_BMU[ie]+=(M_tL_LO_BMU[ie][je]+alphas_muL/4./pi*M_tL_NLO_BMU[ie][je])*CMt_BMU[je];

	return;
}

