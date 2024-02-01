void C_calculator_base2(double C0w[], double C1w[], double mu_W, double C0b[], double C1b[], double mu, struct parameters* param)
/* calculates the LO (C0b) and NLO (C1b) contributions to the Wilson coefficients at scale mu, using the LO (C0w) and NLO (C1w) contributions to the Wilson coefficients at scale mu_W and the parameters of the structure param, in the traditional operator basis */
{
	int ie;
	for(ie=1;ie<=10;ie++) C0b[ie]=C1b[ie]=0.;

	double alphas_muW=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);	
	double eta_mu=alphas_muW/alphas_mu;
	
	double C0w7= C0w[7]-1./3.*C0w[5]-C0w[6]; 
	double C1w7= C1w[7]-1./3.*C1w[5]-C1w[6]; 

	double C0w8= C0w[8]+C0w[5]; 
	double C1w8= C1w[8]+C1w[5]; 

 	C0b[1]= (1./2.*pow(eta_mu,6./23.) -1./2.*pow(eta_mu,-12./23.))*C0w[2];
	C0b[2]= (1./2.*pow(eta_mu,6./23.) +1./2.*pow(eta_mu,-12./23.))*C0w[2];
	C0b[3]= (-1./14.*pow(eta_mu,6./23.) +1./6.*pow(eta_mu,-12./23.) +0.0509*pow(eta_mu,0.4086) -0.1403*pow(eta_mu,-0.4230) -0.01126*pow(eta_mu,-0.8994) +0.0054*pow(eta_mu,0.1456))*C0w[2];
	C0b[4]= (-1./14.*pow(eta_mu,6./23.) -1./6.*pow(eta_mu,-12./23.) +0.0984*pow(eta_mu,0.4086) +0.1214*pow(eta_mu,-0.4230) +0.0156*pow(eta_mu,-0.8994) +0.0026*pow(eta_mu,0.1456))*C0w[2];
	C0b[5]= (-0.0397*pow(eta_mu,0.4086) +0.0117*pow(eta_mu,-0.4230) -0.0025*pow(eta_mu,-0.8994) +0.0304*pow(eta_mu,0.1456))*C0w[2];
	C0b[6]= (0.0335*pow(eta_mu,0.4086) +0.0239*pow(eta_mu,-0.4230) -0.0462*pow(eta_mu,-0.8994) -0.0112*pow(eta_mu,0.1456))*C0w[2];
 
	C0b[7]= pow(eta_mu,16./23.)*C0w[7] + 8./3.*(pow(eta_mu,14./23.)-pow(eta_mu,16./23.))*C0w[8] + C0w[2] * (2.2996*pow(eta_mu,14./23.) -1.0880*pow(eta_mu,16./23.) -3./7.*pow(eta_mu,6./23.) -1./14.*pow(eta_mu,-12./23.) -0.6494*pow(eta_mu,0.4086) -0.0380*pow(eta_mu,-0.4230) -0.0185*pow(eta_mu,-0.8994) -0.0057*pow(eta_mu,0.1456));

	C0b[8]= pow(eta_mu,14./23.)*C0w[8] + C0w[2] * (0.8623*pow(eta_mu,14./23.) -0.9135*pow(eta_mu,0.4086) +0.0873*pow(eta_mu,-0.4230) -0.0571*pow(eta_mu,-0.8994) +0.0209*pow(eta_mu,0.1456));
	
	C0b[9]=C0w[9]+4.*pi/alphas_muW*(-4./33.*(1.-pow(eta_mu,11./23.))+8./87.*(1.-pow(eta_mu,29./23.)))*C0w[2];
	C0b[10]=C0w[10];

	C1b[1]= (C0w[2] *0.8136+1.0197*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)+(C0w[2] *0.7142+2.9524*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.);
	
	C1b[2]= (C0w[2] *0.8136+1.0197*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)-(C0w[2] *0.7142+2.9524*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.);
	
	C1b[3]= (-0.0766*C0w[2]-0.1457*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)+(-0.1455*C0w[2]-0.9841*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.)
	        +(0.1494*eta_mu*C1w[4]-0.8848*C0w[2]+0.2303*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(-0.3726*eta_mu*C1w[4]+0.4137*C0w[2]+1.4672*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(0.0738*eta_mu*C1w[4]-0.0114*C0w[2]+0.0971*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(-0.0173*eta_mu*C1w[4]+0.1722*C0w[2]-0.0213*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	
	C1b[4]= (-0.2353*C0w[2]-0.1457*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)+(-0.0397*C0w[2]+0.9841*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.)
	        +(0.2885*eta_mu*C1w[4]+0.4920*C0w[2]+0.4447*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(0.3224*eta_mu*C1w[4]-0.2758*C0w[2]-1.2696*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(-0.1025*eta_mu*C1w[4]+0.0019*C0w[2]-0.1349*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(-0.0084*eta_mu*C1w[4]-0.1449*C0w[2]-0.0104*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	
	C1b[5]= 0.0397*C0w[2]*pow(eta_mu,6./23.)+0.0926*C0w[2]*pow(eta_mu,-12./23.)
	        +(-0.1163*eta_mu*C1w[4]+0.7342*C0w[2]-0.1792*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(0.0310*eta_mu*C1w[4]-0.1262*C0w[2]-0.1221*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(0.0162*eta_mu*C1w[4]-0.1209*C0w[2]+0.0213*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(-0.0975*eta_mu*C1w[4]-0.1085*C0w[2]-0.1197*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	
	C1b[6]= -0.1191*C0w[2]*pow(eta_mu,6./23.)-0.2778*C0w[2]*pow(eta_mu,-12./23.)
	        +(0.0982*eta_mu*C1w[4]-0.5544*C0w[2]+0.1513*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+(0.0634*eta_mu*C1w[4]+0.1915*C0w[2]-0.2497*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+(0.3026*eta_mu*C1w[4]-0.2744*C0w[2]+0.3983*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+(0.0358*eta_mu*C1w[4]+0.3568*C0w[2]+0.0440*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);

	C1b[7]= pow(eta_mu,39./23.)*C1w7 +8./3.*(pow(eta_mu,37./23.)-pow(eta_mu,39./23.))*C1w8
		+(297664./14283.*pow(eta_mu,16./23.)-7164416./357075.*pow(eta_mu,14./23.)+256868./14283.*pow(eta_mu,37./23.)-6698884./357075.*pow(eta_mu,39./23.))*C0w[8]
		+37208./4761.*(pow(eta_mu,39./23.)-pow(eta_mu,16./23.))*C0w[7]
		+(4661194./816831.*eta_mu*C1w[4]-17.3023*C0w[2]+14.8088*eta_mu*C1w[1]/15.)*pow(eta_mu,14./23.)
		+(-8516./2217.*eta_mu*C1w[4]+8.5027*C0w[2]-10.8090*eta_mu*C1w[1]/15.)*pow(eta_mu,16./23.)
		+ (4.5508*C0w[2]-0.8740*eta_mu*C1w[1]/15.)*pow(eta_mu,6./23.)
		+ (0.7519*C0w[2]+0.4218*eta_mu*C1w[1]/15.)*pow(eta_mu,-12./23.)
		+ (-1.9043*eta_mu*C1w[4]+2.0040*C0w[2]-2.9347*eta_mu*C1w[1]/15.)*pow(eta_mu,0.4086)
		+ (-0.1008*eta_mu*C1w[4]+0.7476*C0w[2]+0.3971*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.4230)
		+ (0.1216*eta_mu*C1w[4]-0.5385*C0w[2]+0.1600*eta_mu*C1w[1]/15.)*pow(eta_mu,-0.8994)
		+ (0.0183*eta_mu*C1w[4]+0.0914*C0w[2]+0.0225*eta_mu*C1w[1]/15.)*pow(eta_mu,0.1456);
	C1b[9]=eta_mu*(C1w[9]+4.*pi/alphas_muW*(-4./33.*(1.-pow(eta_mu,11./23.))+8./87.*(1.-pow(eta_mu,29./23.)))*C1w[2]);
	C1b[10]=eta_mu*C1w[10];
		
	return;	
}

/*----------------------------------------------------------------------*/

void Cprime_calculator(double Cpb[], double complex CQpb[], double mu_W, double mu, struct parameters* param)
{
	int ie;
	for(ie=1;ie<=10;ie++) Cpb[ie]=0.;
	for(ie=1;ie<=2;ie++) CQpb[ie]=0.;
	
	double alphas_muW=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);	
	double eta_mu=alphas_muW/alphas_mu;

	double mass_c_muW=running_mass(param->mass_c,param->mass_c,mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double mass_b_muW=running_mass(param->mass_b,param->mass_b,mu_W,param->mass_top_pole,param->mass_b,param);
	double mass_top_muW=running_mass(param->mtmt,param->mtmt,mu_W,param->mass_top_pole,param->mass_b,param);
	
	/* Wilson coefficients C7 and C8 prime in the SM */ 
	double xt= pow(mass_top_muW/param->mass_W,2.);
	double C7pSM = param->mass_s/mass_b_muW*(-0.5*A0t(xt)-23./36.);
	Cpb[7]=pow(eta_mu,16./23.)*C7pSM;
	double C8pSM = param->mass_s/mass_b_muW*(-0.5*F0t(xt)-1./3.);
	Cpb[8]=pow(eta_mu,14./23.)*C8pSM;
	
	if((param->SM==1)||(param->THDM_model>0)||(param->mass_A02!=0.)||(param->mass_H03!=0.)) return;


 	double sw2=pow(sin(atan(param->gp/param->g2)),2.);
	double MU[4];
		
	MU[1]=param->mass_u;
	MU[2]=mass_c_muW;
	MU[3]=mass_top_muW;

	double epsfac;
	if(param->THDM_model>0) epsfac=1.;
	else epsfac=pow((1.+epsilon_b(param)*param->tan_beta),2.);
	
	double yt= pow(mass_top_muW/param->mass_H,2.);
	double yb= pow(mass_top_muW,2.)/param->mass_c/mass_b_muW;
	double z= pow(param->mass_H/param->mass_W,2.);

	double Gamma_UL[7][4],Gamma_UR[7][4],Gamma_NL[4][4],Gamma_NR[4][4];
	double Gamma_U[7][7], G_aimn[7][3][4][4];
	double X_UL[3][7][4],X_UR[3][7][4],X_NL[3][4][4],X_NR[3][4][4];
	double MD[4],ME[4],VCKM[4][4],Mch[3],MsqU[7],MsqD[7],Msn[4],mintmp;
	double kappa,ag,aY,cosb,sinb,st,ct,alphas_mg;
	double a0a,a0b,a0c,a0Q1,a0Q2,a1,Dp,Dm;
	int ae,be,ce,de,ee,fe,ge,je,ke,me,ne;
	
	VCKM[1][1]=param->Vud;
	VCKM[1][2]=param->Vus;
	VCKM[1][3]=-(param->Vts*param->Vtb+param->Vcs*param->Vcb)/param->Vus; /* Vub from unitarity */
	VCKM[2][1]=param->Vcd;
	VCKM[2][2]=param->Vcs;
	VCKM[2][3]=param->Vcb;
	VCKM[3][1]=param->Vtd;
	VCKM[3][2]=param->Vts;
	VCKM[3][3]=param->Vtb;
	
	sinb=sin(atan(param->tan_beta));
	cosb=cos(atan(param->tan_beta));
	ct=param->stop_mix[2][2];
	st=param->stop_mix[1][2];
	
	MD[1]=param->mass_u;
	MD[2]=param->mass_s;
	MD[3]=mass_b_muW;

	ME[1]=param->mass_e;
	ME[2]=param->mass_mu;
	ME[3]=param->mass_tau;

	Mch[1]=param->mass_cha1;
	Mch[2]=param->mass_cha2;
	
	MsqU[1]=param->mass_upl;
	MsqU[2]=param->mass_chl;
	MsqU[3]=param->mass_t1;
	MsqU[4]=param->mass_upr;
	MsqU[5]=param->mass_chr;
	MsqU[6]=param->mass_t2;
	
	Msn[1]=param->mass_nuel;
	Msn[2]=param->mass_numl;
	Msn[3]=param->mass_nutl;

	
	if(param->THDM_model==0)
	{	alphas_mg=alphas_running(param->mass_gluino,param->mass_top_pole,param->mass_b_pole,param);
		ag=1.-7./12./pi*alphas_mg;
		aY=1.+alphas_mg/4./pi;
		
		kappa=1./(param->g2*param->g2*param->Vtb*param->Vts);
					if(param->sU_mix[1][1]*param->sU_mix[1][2]*param->sU_mix[1][3]*param->sU_mix[1][4]*param->sU_mix[1][5]*param->sU_mix[1][6]!=0.)
		{
		
			for(je=1;je<=5;je++)
			{
				for(ie=je+1;ie<=6;ie++) if (MsqU[je]>MsqU[ie])
				{
					mintmp=MsqU[je];
					MsqU[je]=MsqU[ie];
					MsqU[ie]=mintmp;
				}
			}
		
		
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UL[ae][ie]=param->sU_mix[ae][ie];
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UR[ae][ie]=param->sU_mix[ae][ie+3];
		}
		else
		{
			Gamma_UL[1][1]=1.;
			Gamma_UL[2][1]=0.;
			Gamma_UL[3][1]=0.;
			Gamma_UL[4][1]=0.;
			Gamma_UL[5][1]=0.;
			Gamma_UL[6][1]=0.;
			Gamma_UL[1][2]=0.;
			Gamma_UL[2][2]=1.;
			Gamma_UL[3][2]=0.;
			Gamma_UL[4][2]=0.;
			Gamma_UL[5][2]=0.;
			Gamma_UL[6][2]=0.;
			Gamma_UL[1][3]=0.;
			Gamma_UL[2][3]=0.;
			Gamma_UL[3][3]=ct;
			Gamma_UL[4][3]=0.;
			Gamma_UL[5][3]=0.;
			Gamma_UL[6][3]=-st;
		
			Gamma_UR[1][1]=0.;
			Gamma_UR[2][1]=0.;
			Gamma_UR[3][1]=0.;
			Gamma_UR[4][1]=1.;
			Gamma_UR[5][1]=0.;
			Gamma_UR[6][1]=0.;
			Gamma_UR[1][2]=0.;
			Gamma_UR[2][2]=0.;
			Gamma_UR[3][2]=0.;
			Gamma_UR[4][2]=0.;
			Gamma_UR[5][2]=1.;
			Gamma_UR[6][2]=0.;
			Gamma_UR[1][3]=0.;
			Gamma_UR[2][3]=0.;
			Gamma_UR[3][3]=st;
			Gamma_UR[4][3]=0.;
			Gamma_UR[5][3]=0.;
			Gamma_UR[6][3]=ct;
		}

		for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++)
		{
			Gamma_U[ae][ie]=Gamma_UL[ae][ie];
			Gamma_U[ae][ie+3]=Gamma_UR[ae][ie];
		}
						
		for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++) if(ae==be) Gamma_NL[ae][be]=Gamma_NR[ae][be]=1.; else Gamma_NL[ae][be]=Gamma_NR[ae][be]=0.;
				
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
		{
			X_UL[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UL[ie][ae][be]+=-param->g2*(ag*param->charg_Vmix[ie][1]*Gamma_UL[ae][ce]
			-aY*param->charg_Vmix[ie][2]*Gamma_UR[ae][ce]*MU[ce]/(sqrt(2.)*param->mass_W*sinb))*VCKM[ce][be];
		
			X_UR[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UR[ie][ae][be]+=param->g2*aY*param->charg_Umix[ie][2]*Gamma_UL[ae][ce]*VCKM[ce][be]*MD[be]/(sqrt(2)*param->mass_W*cosb);
		}
	
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++)
		{
			X_NL[ie][ae][be]=-param->g2*param->charg_Vmix[ie][1]*Gamma_NL[ae][be];
			X_NR[ie][ae][be]=param->g2*param->charg_Umix[ie][2]*Gamma_NL[ae][be]*ME[be]/(sqrt(2.)*param->mass_W*cosb);
		}
		
		for(ae=1;ae<=6;ae++) for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
		{
			G_aimn[ae][ie][me][ne]=0.5/sqrt(2.)*(sqrt(2.)*param->mass_W*param->charg_Vmix[ie][1]*Gamma_UL[ae][ne]*ag-MU[ne]*param->charg_Vmix[ie][2]*Gamma_UR[ae][ne]*aY)*(VCKM[me][3]*VCKM[ne][2]/VCKM[3][3]/VCKM[3][2]);
		}

	}
	
	/* Wilson coefficients C7 and C8 prime */ 
	double ld=-param->tan_beta;
	double C7pH=param->mass_s*mass_b_muW/mass_top_muW/mass_top_muW*1./3.*ld*ld*F7_1(yt);
	double C8pH=param->mass_s*mass_b_muW/mass_top_muW/mass_top_muW*1./3.*ld*ld*F8_1(yt);

	double C7pcharg=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C7pcharg+=pow(param->mass_W/Mch[ie],2.)*(X_UR[ie][ae][2]*X_UR[ie][ae][3]*h10(pow(MsqU[ae]/Mch[ie],2.)) + Mch[ie]/mass_b_muW*X_UR[ie][ae][2]*X_UL[ie][ae][3]*h20(pow(MsqU[ae]/Mch[ie],2.)));		
	C7pcharg*=-0.5*kappa; 
		
	double C8pcharg=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) C8pcharg+=pow(param->mass_W/Mch[ie],2.)*(X_UR[ie][ae][2]*X_UR[ie][ae][3]*h50(pow(MsqU[ae]/Mch[ie],2.)) + Mch[ie]/mass_b_muW*X_UR[ie][ae][2]*X_UL[ie][ae][3]*h60(pow(MsqU[ae]/Mch[ie],2.)));		
	C8pcharg*=-0.5*kappa; 

	Cpb[7]+=pow(eta_mu,16./23.)*(C7pH+C7pcharg);
	Cpb[8]+=pow(eta_mu,14./23.)*(C8pH+C8pcharg);
	/* printf("C7p=%.5e\n",Cpb[7]); */
	/* printf("C8p=%.5e\n",Cpb[8]); */
	
	
	/* Wilson coefficients C9 and C10 prime */ 	
	double C10pH = -mass_b_muW*param->mass_s*(param->tan_beta*param->tan_beta/8./param->mass_W/param->mass_W
	+pow(param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_H,2.))*f20(yt)/sw2;
	
	double C9pH =(4.*sw2-1.)*C10pH - param->mass_s*mass_b_muW/mass_top_muW/mass_top_muW*D9H0(yt,ld);
	
	
	double B10pc=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++) B10pc+=-X_UR[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(0.5*X_NR[ie][be][2]*X_NR[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)) +X_NL[ie][be][2]*X_NL[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
	B10pc*=kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;
		
	double C9pc=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) C9pc+=X_UR[je][ae][2]*X_UR[ie][ae][3]*(2.*fabs(Mch[je]/Mch[ie])*f30(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))*param->charg_Vmix[je][1]*param->charg_Vmix[ie][1] -f40(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.))*param->charg_Umix[je][1]*param->charg_Umix[ie][1]);
		
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=3;ce++) C9pc+=-X_UR[ie][be][2]*X_UR[ie][ae][3]*f40(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))*Gamma_UR[be][ce]*Gamma_UR[ae][ce];
	C9pc*=-kappa/8.;
		
	double C10pcharg=(B10pc-C9pc)/sw2;
	
	double B9pc=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++) B9pc+=X_UR[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(0.5*X_NR[ie][be][2]*X_NR[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)) -X_NL[ie][be][2]*X_NL[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
	B9pc*=kappa*param->mass_W*param->mass_W/2./param->g2/param->g2;
	
	double D9pc=0.;
	for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) D9pc+=pow(param->mass_W/Mch[ie],2.)*X_UR[ie][ae][2]*X_UR[ie][ae][3]*h30(pow(MsqU[ae]/Mch[ie],2.));
	D9pc*=kappa;

	double C9pcharg=(1.-4.*sw2)/sw2*C9pc-B9pc/sw2-D9pc;
	
	Cpb[9]+=C9pH+C9pcharg;
	Cpb[10]+=C10pH+C10pcharg;	
	/* printf("C9p=%.5e\n",Cpb[9]); */
	/* printf("C10p=%.5e\n",Cpb[10]); */


	/* Wilson coefficients CQ1 and CQ2 prime */ 
	double NQ1pH=-param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*xt*f30(xt,z);
	
	double BQ1pH=param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*f70(xt,z);
	
	double complex CQ1pH=(NQ1pH+BQ1pH)*param->mass_s/sw2;
	
	double complex CQ2pH=CQ1pH;
	
	
	double BQ1pc1=0.;
	double BQ1pc2=0.;
	double BQ2pc1=0.;
	double BQ2pc2=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
	{ 	BQ1pc1+=X_UR[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NL[ie][be][2]*X_NR[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))
	); 	BQ1pc2+=X_UR[je][ae][2]*X_UL[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NR[ie][be][2]*X_NL[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
	}
	
	double BQ1pc=(BQ1pc1+BQ1pc2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	double BQ2pc=(BQ1pc1-BQ1pc2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;

	
	double NQ1pc=0.;
	double NQ2pc=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++)
		{ Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]+Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]-Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
		}
		a0a=-(fabs(Mch[ie]/Mch[je])*f30(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))*param->charg_Umix[ie][2]*param->charg_Vmix[je][1])*kron(ae,be);

a0b=-(f40(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))*param->charg_Umix[je][2]*param->charg_Vmix[ie][1])*kron(ae,be);
		a0c=1./param->mass_W*f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))*kron(ie,je);
		a0Q1=a0a+a0b+Dp*a0c;
		a0Q2=-a0a+a0b+Dm*a0c;
		a1=Mch[ie]/sqrt(2.)/param->mass_W*f80(pow(MsqU[ae]/Mch[ie],2.))*kron(ie,je)*kron(ae,be);
		
		NQ1pc+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q1+a1*param->tan_beta);
		NQ2pc+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q2+a1*param->tan_beta);

	}
	NQ1pc*=param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*param->mass_s/sw2;	NQ2pc*=param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*param->mass_s/sw2;
	
	double complex CQ1pcharg=NQ1pc+BQ1pc;
	CQpb[1]=CQ1pH+CQ1pcharg;
	CQpb[1]/=epsfac;
	

	double complex CQ2pcharg=NQ2pc+BQ2pc;
	CQpb[2]=CQ2pH+CQ2pcharg;
	CQpb[2]/=epsfac;

	int nf=5;
	double beta0 = 11.-2./3.*nf;
	CQpb[1]*=pow(eta_mu,-4./beta0);
	CQpb[2]*=pow(eta_mu,-4./beta0);
	/* printf("CQ1p=%.5e\n",CQpb[1]);
	printf("CQ2p=%.5e\n",CQpb[2]); */

	return;
}


/*-----------------------------------------------------------------------------*/

	void CQ_calculator(double complex CQ0b[], double complex CQ1b[], double mu_W, double mu, struct parameters* param)
{
	int ie;
	for(ie=1;ie<=2;ie++) CQ0b[ie]=CQ1b[ie]=0.;
	
	if(param->SM==1) return;

 	double sw2=pow(sin(atan(param->gp/param->g2)),2.);
	double MU[4];
	
	double mass_c_muW=running_mass(param->mass_c,param->mass_c,mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double mass_b_muW=running_mass(param->mass_b,param->mass_b,mu_W,param->mass_top_pole,param->mass_b,param);
	double mass_top_muW=running_mass(param->mtmt,param->mtmt,mu_W,param->mass_top_pole,param->mass_b,param);
	
	MU[1]=param->mass_u;
	MU[2]=mass_c_muW;
	MU[3]=mass_top_muW;

	double alphas_muW=alphas_running(mu_W,param->mass_top_pole,param->mass_b_pole,param);
	double alphas_mu=alphas_running(mu,param->mass_top_pole,param->mass_b_pole,param);	
	double eta_mu=alphas_muW/alphas_mu;

	int nf=5;
	double beta0 = 11.-2./3.*nf;
	
	double complex CQ1H_0,CQ2H_0;

	if(param->THDM_model==0) 
	{	
		param->lambda_u[3][3]=1./param->tan_beta;
		param->lambda_d[2][2]=-param->tan_beta;
		param->lambda_d[3][3]=-param->tan_beta;
		param->lambda_l[2][2]=-param->tan_beta;
	}

	/* Wilson coefficients CQ1 et CQ2 in 2HDM */ 
 
	double complex CQ1box=-param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]*pow(1./2./param->mass_W,2.)/sw2*Bplus(param->mass_H*param->mass_H/param->mass_W/param->mass_W,MU[3]*MU[3]/param->mass_W/param->mass_W)/4.;

	double complex CQ2box=-CQ1box;
		
	double Pplus=-D2(param->mass_H*param->mass_H/param->mass_W/param->mass_W,MU[3]*MU[3]/param->mass_W/param->mass_W)*MU[3]*MU[3]/param->mass_W/param->mass_W;
	
	double complex CQ1peng1=-param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]/4./sw2*Pplus*(sin(param->alpha)*sin(param->alpha)/param->mass_h0/param->mass_h0+cos(param->alpha)*cos(param->alpha)/param->mass_H0/param->mass_H0)/4.;
	
	double complex CQ2peng1=param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]/4./sw2*Pplus/param->mass_A0/param->mass_A0/4.;
		
	double complex CQ1peng2=param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]/4./sw2*Pplus*(sin(param->alpha)*sin(param->alpha)/param->mass_h0/param->mass_h0*(param->mass_H*param->mass_H-param->mass_h0*param->mass_h0)/param->mass_W/param->mass_W+cos(param->alpha)*cos(param->alpha)/param->mass_H0/param->mass_H0*(param->mass_H*param->mass_H-param->mass_H0*param->mass_H0)/param->mass_W/param->mass_W)/4.;
	
	double complex CQ2peng2=-param->mass_mu*(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3])*param->lambda_l[2][2]/4./sw2*Pplus/param->mass_A0/param->mass_A0*(param->mass_H*param->mass_H-param->mass_A0*param->mass_A0)/param->mass_W/param->mass_W/4.;
 
	
	double complex CQ1self=-param->mass_mu*param->lambda_d[3][3]*param->lambda_l[2][2]/4./sw2*(param->mass_H*param->mass_H/param->mass_W/param->mass_W+((param->lambda_d[3][3]+param->mass_s/mass_b_muW*param->lambda_d[2][2])*param->lambda_u[3][3]))*Pplus*(sin(param->alpha)*sin(param->alpha)/param->mass_h0/param->mass_h0+cos(param->alpha)*cos(param->alpha)/param->mass_H0/param->mass_H0)/4.;
	
	double complex CQ2self=param->mass_mu*param->lambda_d[3][3]*param->lambda_l[2][2]/4./sw2*(param->mass_H*param->mass_H/param->mass_W/param->mass_W+((param->lambda_d[3][3]-param->mass_s/mass_b_muW*param->lambda_d[2][2])*param->lambda_u[3][3]))*Pplus/param->mass_A0/param->mass_A0/4.;
		
	
	CQ1box *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/param->mass_s*param->lambda_u[3][3]);
	
	CQ2box *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);
	
	CQ1peng1 *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);
	
	CQ2peng1 *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);

	CQ1peng2 *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);
	
	CQ2peng2 *= 1.+param->mass_s/mass_b_muW*(param->lambda_d[2][2]-MU[3]/param->mass_s*param->lambda_u[3][3])/(param->lambda_d[3][3]-MU[3]/mass_b_muW*param->lambda_u[3][3]);
	
	CQ1self *= 1.+param->mass_s/mass_b_muW*param->lambda_d[2][2]/param->lambda_d[3][3];

	CQ2self *= 1.+param->mass_s/mass_b_muW*param->lambda_d[2][2]/param->lambda_d[3][3];
	
	CQ1H_0=(CQ1box+CQ1peng1+CQ1peng2+CQ1self);
	CQ2H_0=(CQ2box+CQ2peng1+CQ2peng2+CQ2self);
		
	if(param->THDM_model>0)
	{
		CQ0b[1]=(CQ1H_0)*mass_b_muW/sw2;
		CQ0b[2]=(CQ2H_0)*mass_b_muW/sw2;
	
		CQ0b[1]*=pow(eta_mu,-4./beta0);
		CQ0b[2]*=pow(eta_mu,-4./beta0);

		return;
	}
	
	double epsfac=pow((1.+epsilon_b(param)*param->tan_beta),2.);
	
	double xt= pow(mass_top_muW/param->mass_W,2.);
	double yt= pow(mass_top_muW/param->mass_H,2.);
	double yb= pow(mass_top_muW,2.)/param->mass_c/mass_b_muW;
	double z= pow(param->mass_H/param->mass_W,2.);

	double Gamma_UL[7][4],Gamma_UR[7][4],Gamma_NL[4][4],Gamma_NR[4][4];
	double Gamma_U[7][7], G_aimn[7][3][4][4],I_LR[7][7],P_U[7][7];
	double X_UL[3][7][4],X_UR[3][7][4],X_NL[3][4][4],X_NR[3][4][4];
	double MD[4],ME[4],VCKM[4][4],Mch[3],MsqU[7],MsqD[7],Msn[4],mintmp;
	double kappa,ag,aY,cosb,sinb,st,ct,alphas_mg,temp;
	double a0,a0Q1,a0Q2,a1,a0a,a0b,a0c,a0p,a2p,Dp,Dm,Dpm;
	int ae,be,ce,de,ee,fe,ge,je,ke,me,ne;

		
	VCKM[1][1]=param->Vud;
	VCKM[1][2]=param->Vus;
	VCKM[1][3]=-(param->Vts*param->Vtb+param->Vcs*param->Vcb)/param->Vus; /* Vub from unitarity */
	VCKM[2][1]=param->Vcd;
	VCKM[2][2]=param->Vcs;
	VCKM[2][3]=param->Vcb;
	VCKM[3][1]=param->Vtd;
	VCKM[3][2]=param->Vts;
	VCKM[3][3]=param->Vtb;
		
	sinb=sin(atan(param->tan_beta));
	cosb=cos(atan(param->tan_beta));
	ct=param->stop_mix[2][2];
	st=param->stop_mix[1][2];
		
	MD[1]=param->mass_u;
	MD[2]=param->mass_s;
	MD[3]=mass_b_muW;

	ME[1]=param->mass_e;
	ME[2]=param->mass_mu;
	ME[3]=param->mass_tau;

	Mch[1]=param->mass_cha1;
	Mch[2]=param->mass_cha2;
		
	MsqU[1]=param->mass_upl;
	MsqU[2]=param->mass_chl;
	MsqU[3]=param->mass_t1;
	MsqU[4]=param->mass_upr;
	MsqU[5]=param->mass_chr;
	MsqU[6]=param->mass_t2;
		
	Msn[1]=param->mass_nuel;
	Msn[2]=param->mass_numl;
	Msn[3]=param->mass_nutl;

	if((param->mass_A02==0.)&&(param->mass_H03==0.))
	{
		alphas_mg=alphas_running(param->mass_gluino,param->mass_top_pole,param->mass_b_pole,param);
		ag=1.-7./12./pi*alphas_mg;
		aY=1.+alphas_mg/4./pi;
	
		kappa=1./(param->g2*param->g2*param->Vtb*param->Vts);
		
	if(param->sU_mix[1][1]*param->sU_mix[1][2]*param->sU_mix[1][3]*param->sU_mix[1][4]*param->sU_mix[1][5]*param->sU_mix[1][6]!=0.)
		{
		
			for(je=1;je<=5;je++)
			{
				for(ie=je+1;ie<=6;ie++) if (MsqU[je]>MsqU[ie])
				{
					mintmp=MsqU[je];
					MsqU[je]=MsqU[ie];
					MsqU[ie]=mintmp;
				}
			}
		
		
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UL[ae][ie]=param->sU_mix[ae][ie];
			for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++) Gamma_UR[ae][ie]=param->sU_mix[ae][ie+3];
		}
		else
		{
			Gamma_UL[1][1]=1.;
			Gamma_UL[2][1]=0.;
			Gamma_UL[3][1]=0.;
			Gamma_UL[4][1]=0.;
			Gamma_UL[5][1]=0.;
			Gamma_UL[6][1]=0.;
			Gamma_UL[1][2]=0.;
			Gamma_UL[2][2]=1.;
			Gamma_UL[3][2]=0.;
			Gamma_UL[4][2]=0.;
			Gamma_UL[5][2]=0.;
			Gamma_UL[6][2]=0.;
			Gamma_UL[1][3]=0.;
			Gamma_UL[2][3]=0.;
			Gamma_UL[3][3]=ct;
			Gamma_UL[4][3]=0.;
			Gamma_UL[5][3]=0.;
			Gamma_UL[6][3]=-st;
		
			Gamma_UR[1][1]=0.;
			Gamma_UR[2][1]=0.;
			Gamma_UR[3][1]=0.;
			Gamma_UR[4][1]=1.;
			Gamma_UR[5][1]=0.;
			Gamma_UR[6][1]=0.;
			Gamma_UR[1][2]=0.;
			Gamma_UR[2][2]=0.;
			Gamma_UR[3][2]=0.;
			Gamma_UR[4][2]=0.;
			Gamma_UR[5][2]=1.;
			Gamma_UR[6][2]=0.;
			Gamma_UR[1][3]=0.;
			Gamma_UR[2][3]=0.;
			Gamma_UR[3][3]=st;
			Gamma_UR[4][3]=0.;
			Gamma_UR[5][3]=0.;
			Gamma_UR[6][3]=ct;
		}

		for(ae=1;ae<=6;ae++) for(ie=1;ie<=3;ie++)
		{
			Gamma_U[ae][ie]=Gamma_UL[ae][ie];
			Gamma_U[ae][ie+3]=Gamma_UR[ae][ie];
		}
		
		for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) I_LR[ae][be]=0.;
		for(ae=1;ae<=3;ae++) I_LR[ae][ae]=1.;
		for(ae=4;ae<=6;ae++) I_LR[ae][ae]=-1.;
		
		for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) for(de=1;de<=6;de++) P_U[ae][be]=Gamma_U[ae][ce]*I_LR[ce][de]*Gamma_U[be][de];
		
		for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++) if(ae==be) Gamma_NL[ae][be]=Gamma_NR[ae][be]=1.; else Gamma_NL[ae][be]=Gamma_NR[ae][be]=0.;
				
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
		{
			X_UL[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UL[ie][ae][be]+=-param->g2*(ag*param->charg_Vmix[ie][1]*Gamma_UL[ae][ce]-aY*param->charg_Vmix[ie][2]*Gamma_UR[ae][ce]*MU[ce]/(sqrt(2.)*param->mass_W*sinb))*VCKM[ce][be];
		
			X_UR[ie][ae][be]=0.;
			for(ce=1;ce<=3;ce++) X_UR[ie][ae][be]+=param->g2*aY*param->charg_Umix[ie][2]*Gamma_UL[ae][ce]*VCKM[ce][be]*MD[be]/(sqrt(2)*param->mass_W*cosb);
		}
	
		for(ie=1;ie<=2;ie++) for(ae=1;ae<=3;ae++) for(be=1;be<=3;be++)
		{
			X_NL[ie][ae][be]=-param->g2*param->charg_Vmix[ie][1]*Gamma_NL[ae][be];
			X_NR[ie][ae][be]=param->g2*param->charg_Umix[ie][2]*Gamma_NL[ae][be]*ME[be]/(sqrt(2.)*param->mass_W*cosb);
		}
		
		for(ae=1;ae<=6;ae++) for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
		{
			G_aimn[ae][ie][me][ne]=0.5/sqrt(2.)*(sqrt(2.)*param->mass_W*param->charg_Vmix[ie][1]*Gamma_UL[ae][ne]*ag-MU[ne]*param->charg_Vmix[ie][2]*Gamma_UR[ae][ne]*aY)*(VCKM[me][3]*VCKM[ne][2]/VCKM[3][3]/VCKM[3][2]);
		}
	
	/* LO */
	
	double NQ10H=-param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*xt*f30(xt,z);
	double BQ10H=param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*f70(xt,z);
	CQ1H_0=(NQ10H+BQ10H)*mass_b_muW/sw2; 
	CQ2H_0=-CQ1H_0;
	
	
	double BQ10c1=0.;
	double BQ10c2=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++)
	{ BQ10c1+=X_UL[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NR[ie][be][2]*X_NL[je][be][2]*f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
BQ10c2+=X_UL[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NL[ie][be][2]*X_NR[je][be][2]*fabs(Mch[je]/Mch[ie])*f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.)));
	}
	double BQ10c=(BQ10c1+BQ10c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	double BQ20c=-(BQ10c1-BQ10c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	
	
	double NQ10c=0.;
	double NQ20c=0.;
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++) 
		{
			Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]+Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
			Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]-Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
		}
			a0a=-(fabs(Mch[ie]/Mch[je])*f30(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))*param->charg_Umix[ie][2]*param->charg_Vmix[je][1])*kron(ae,be);
a0b=-(f40(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))*param->charg_Umix[je][2]*param->charg_Vmix[ie][1])*kron(ae,be);
		a0c=1./param->mass_W*f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))*kron(ie,je);
		a0Q1=a0a+a0b+Dp*a0c;
		a0Q2=-a0a+a0b+Dm*a0c;
		
		a1=Mch[ie]/sqrt(2.)/param->mass_W*f80(pow(MsqU[ae]/Mch[ie],2.))*kron(ie,je)*kron(ae,be);
		NQ10c+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q1+a1*param->tan_beta);
		NQ20c+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q2+a1*param->tan_beta);
	}
	NQ10c*=param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;
        NQ20c*=-param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;


	double complex CQ1charg_0=NQ10c+BQ10c;
	double complex CQ2charg_0=NQ20c+BQ20c;
	
	CQ0b[1]=CQ1H_0+CQ1charg_0;
	CQ0b[1]/=epsfac;
		
	CQ0b[2]=CQ2H_0+CQ2charg_0;
	CQ0b[2]/=epsfac;
	
	/* NLO - Charged Higgs */
	
	double NQ11H=-param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*(f141(xt,z)+8.*xt*(f30(xt,z)+xt*(f30(xt*1.0001,z)-f30(xt*0.9999,z))/0.0002)*log(mu_W*mu_W/mass_top_muW/mass_top_muW));
	double BQ11H=param->mass_mu*param->tan_beta*param->tan_beta/4./param->mass_W/param->mass_W*(f111(xt,z)+8.*(f70(xt*1.0001,z)-f70(xt*0.9999,z))/0.0002*log(mu_W*mu_W/mass_top_muW/mass_top_muW));
	double CQ1H_1=(NQ11H+BQ11H)*mass_b_muW/sw2;
	double CQ2H_1=-CQ1H_1;

	/* NLO - charginos */
	
	double BQ11c1=0.;
	double BQ11c2=0.;
	
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++) BQ11c1+=X_UL[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NR[ie][be][2]*X_NL[je][be][2]*(f121(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))+4.*(f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*1.0001,pow(Msn[be]/Mch[ie],2.))-f50(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*0.9999,pow(Msn[be]/Mch[ie],2.)))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.))));
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=3;be++) BQ11c2+=X_UL[je][ae][2]*X_UR[ie][ae][3]/Mch[ie]/Mch[ie]*(X_NL[ie][be][2]*X_NR[je][be][2]*fabs(Mch[je]/Mch[ie])*(f131(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(Msn[be]/Mch[ie],2.))+4.*(f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*1.0001,pow(Msn[be]/Mch[ie],2.))-f60(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.)*0.9999,pow(Msn[be]/Mch[ie],2.)))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.))));
	
	double BQ11c=(BQ11c1+BQ11c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	double BQ21c=-(BQ11c1-BQ11c2)*kappa*param->mass_W*param->mass_W/2./param->g2/param->g2/sw2;
	
	
	double NQ11c=0.;
	double NQ21c=0.;

	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++)
		{ 	Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]+Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
	Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[be][fe]-Gamma_UL[ae][fe]*Gamma_UR[be][fe]);
		}
			a0a=-(fabs(Mch[ie]/Mch[je])*(f181(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))+4.*(f30(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.)*1.0001)-f30(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.)*0.9999))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.)))*param->charg_Umix[ie][2]*param->charg_Vmix[je][1])*kron(ae,be);
					a0b=-((f191(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.))+4.*(f40(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.)*1.0001)-f40(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.)*0.9999))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.)))*param->charg_Umix[je][2]*param->charg_Vmix[ie][1])*kron(ae,be);
					a0c=1./param->mass_W*(f171(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))+4.*(f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.))+(f30(pow(MsqU[ae]/Mch[ie],2.)*1.0001,pow(MsqU[be]/Mch[ie],2.))-f30(pow(MsqU[ae]/Mch[ie],2.)*0.9999,pow(MsqU[be]/Mch[ie],2.)))/0.0002+(f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.)*1.0001)-f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.)*0.9999))/0.0002)*log(pow(mass_top_muW/MsqU[ae],2.)))*kron(ie,je);
	
		a0Q1=a0a+a0b+Dp*a0c;
		a0Q2=-a0a+a0b+Dm*a0c;
			a0p=4.*G_aimn[ae][ie][me][ne]/param->mass_W/(VCKM[me][3]*VCKM[ne][2]/VCKM[3][3]/VCKM[3][2])/param->charg_Umix[je][2]*f151(pow(MsqU[ae]/Mch[ie],2.))*kron(ie,je)*kron(ae,be)*kron(me,ne);
				a1=Mch[ie]/sqrt(2.)/param->mass_W*(f161(pow(MsqU[ae]/Mch[ie],2.))+4.*(f80(pow(MsqU[ae]/Mch[ie],2.)*1.0001)-f80(pow(MsqU[ae]/Mch[ie],2.)*0.9999))/0.0002*log(pow(mass_top_muW/MsqU[ae],2.)))*kron(ie,je)*kron(ae,be);
		a2p=Gamma_UL[be][me]*(VCKM[me][3]*VCKM[ne][2]/VCKM[3][3]/VCKM[3][2])*param->charg_Umix[je][2]/2./param->mass_W*f151(pow(MsqU[ae]/Mch[ie],2.))*kron(ie,je)*kron(ae,be)*kron(me,ne);
		
		NQ11c+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q1+a1*param->tan_beta)
		+G_aimn[ae][ie][me][ne]*param->charg_Umix[je][2]*a0p
		+Gamma_UL[be][me]*param->charg_Umix[je][2]*a2p*pow(param->mass_s*param->tan_beta,2.);	
		NQ21c+=G_aimn[ae][ie][me][ne]*Gamma_UL[be][me]*param->charg_Umix[je][2]*(a0Q2+a1*param->tan_beta)
		+G_aimn[ae][ie][me][ne]*param->charg_Umix[je][2]*a0p
		+Gamma_UL[be][me]*param->charg_Umix[je][2]*a2p*pow(param->mass_s*param->tan_beta,2.);
	}
		NQ11c*=param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;
	NQ21c*=-param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;
	
	double complex CQ1charg_1=NQ11c+BQ11c;
	
	if(cabs(CQ1charg_1)*alphas_mu/4./pi>cabs(CQ0b[1])) CQ1charg_1*=cabs(CQ0b[1])/cabs(CQ1charg_1)*4.*pi/alphas_mu;
	if(cabs(CQ1H_1)*alphas_mu/4./pi>cabs(CQ0b[1])) CQ1H_1*=cabs(CQ0b[1])/cabs(CQ1H_1)*4.*pi/alphas_mu;
	CQ1b[1]=CQ1H_1+CQ1charg_1;
	
	CQ1b[1]/=epsfac;
	
	double complex CQ2charg_1=NQ21c+BQ21c;
	
	if(cabs(CQ2charg_1)*alphas_mu/4./pi>cabs(CQ0b[2])) CQ2charg_1*=cabs(CQ0b[2])/cabs(CQ2charg_1)*4.*pi/alphas_mu;
	if(cabs(CQ2H_1)*alphas_mu/4./pi>cabs(CQ0b[2])) CQ2H_1*=cabs(CQ0b[2])/cabs(CQ2H_1)*4.*pi/alphas_mu;
	CQ1b[2]=CQ2H_1+CQ2charg_1;
		
	CQ1b[2]/=epsfac;
	
	
	/* Wilson coefficient CQ1 */ 
	/* NLO  - four points */
	double BQ11f1=0.;
	double BQ11f2=0.;
	
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(fe=1;fe<=3;fe++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++)
	{ BQ11f1+=-X_UL[je][be][2]*X_UR[ie][ae][3]*pow(param->mass_W/Mch[ie],2.)*P_U[ae][ce]*MsqU[ce]/Mch[ie]*P_U[ce][be]*(1.+log(pow(mass_top_muW/MsqU[ce],2.)))	*(f90(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.),pow(Msn[fe]/Mch[ie],2.))*X_NR[ie][fe][2]*X_NL[je][fe][2]);
		BQ11f2+=-X_UL[je][be][2]*X_UR[ie][ae][3]*pow(param->mass_W/Mch[ie],2.)*P_U[ae][ce]*MsqU[ce]/Mch[ie]*P_U[ce][be]*(1.+log(pow(mass_top_muW/MsqU[ce],2.)))	*(fabs(Mch[je]/Mch[ie])*f100(pow(Mch[je]/Mch[ie],2.),pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.),pow(Msn[fe]/Mch[ie],2.))*X_NL[ie][fe][2]*X_NR[je][fe][2]);
	}
		
	double BQ11f=(BQ11f1+BQ11f2)*2./3.*kappa/param->g2/param->g2/sw2;
	double BQ21f=-(BQ11f1-BQ11f2)*2./3.*kappa/param->g2/param->g2/sw2;
	
	double NQ11f=0.;
	double NQ21f=0.;

	for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(de=1;de<=6;de++) for(ke=1;ke<=6;ke++) NQ11f+=G_aimn[ae][ie][me][ne]*Gamma_UL[de][me]*param->charg_Umix[ie][2]*P_U[ae][ke]*MsqU[ke]/Mch[ie]*P_U[ke][de]*(1.+log(pow(mass_top_muW/MsqU[ke],2.)))*param->tan_beta*Mch[ie]/sqrt(2.)*f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[de]/Mch[ie],2.));

	NQ21f=NQ11f;

	for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(be=1;be<=6;be++) for(ce=1;ce<=6;ce++) for(ke=1;ke<=6;ke++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++) 
		{		Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[be][fe]*Gamma_UL[ce][fe]+Gamma_UL[be][fe]*Gamma_UR[ce][fe]); 
Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[be][fe]*Gamma_UL[ce][fe]-Gamma_UL[be][fe]*Gamma_UR[ce][fe]); 
		}
		temp=G_aimn[ae][ie][me][ne]*Gamma_UL[ce][me]*param->charg_Umix[ie][2]*P_U[ae][ke]*MsqU[ke]/Mch[ie]*P_U[ke][be]*(1.+log(pow(mass_top_muW/MsqU[ke],2.)))*f60(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[be]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.));	
		NQ11f+=Dp*temp;
		NQ21f+=Dm*temp;
	}
	
	for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(ce=1;ce<=6;ce++)  for(de=1;de<=6;de++) for(ke=1;ke<=6;ke++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++)
		{ Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[ce][fe]+Gamma_UL[ae][fe]*Gamma_UR[ce][fe]); Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ae][fe]*Gamma_UL[ce][fe]-Gamma_UL[ae][fe]*Gamma_UR[ce][fe]); 
		}
	
	temp=G_aimn[ae][ie][me][ne]*Gamma_UL[de][me]*param->charg_Umix[ie][2]*P_U[ce][ke]*MsqU[ke]/Mch[ie]*P_U[ke][de]*(1.+log(pow(mass_top_muW/MsqU[ke],2.)))*f60(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[ce]/Mch[ie],2.),pow(MsqU[de]/Mch[ie],2.));	
	NQ11f+=temp*Dp;
	NQ21f+=temp*Dm;
	}
	
	for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(de=1;de<=6;de++) for(ke=1;ke<=6;ke++)
	{
		temp=-G_aimn[ae][ie][me][ne]*Gamma_UL[de][me]*param->charg_Umix[je][2]*P_U[ae][ke]*MsqU[ke]/Mch[je]*P_U[ke][de]*(1.+log(pow(mass_top_muW/MsqU[ke],2.)))*param->mass_W*(fabs(Mch[ie]/Mch[je])*f60(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.),pow(MsqU[de]/Mch[je],2.))*param->charg_Umix[ie][2]*param->charg_Vmix[je][1]+f50(pow(Mch[ie]/Mch[je],2.),pow(MsqU[ae]/Mch[je],2.),pow(MsqU[de]/Mch[je],2.))*param->charg_Umix[je][2]*param->charg_Vmix[ie][1]); 
	NQ11f+=temp;
	NQ21f+=-temp;
	}

	
	for(ie=1;ie<=2;ie++) for(me=1;me<=3;me++) for(ne=1;ne<=3;ne++) for(ae=1;ae<=6;ae++) for(de=1;de<=6;de++)for(ee=1;ee<=6;ee++) for(ge=1;ge<=6;ge++)
	{
		Dp=0.;
		Dm=0.;
		for(fe=1;fe<=3;fe++)
		{ Dp+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ee][fe]*Gamma_UL[ge][fe]+Gamma_UL[ee][fe]*Gamma_UR[ge][fe]); 
Dm+=MU[fe]/sqrt(2.)/Mch[ie]*param->mu_Q*(Gamma_UR[ee][fe]*Gamma_UL[ge][fe]-Gamma_UL[ee][fe]*Gamma_UR[ge][fe]); 
		}
	
	temp=	-G_aimn[ae][ie][me][ne]*Gamma_UL[de][me]*param->charg_Umix[ie][2]*P_U[ae][ee]*(1.+log(pow(mass_top_muW/MsqU[ge],2.))-f110(pow(MsqU[ee]/Mch[ie],2.),pow(MsqU[ge]/Mch[ie],2.)))*P_U[ge][de]*f30(pow(MsqU[ae]/Mch[ie],2.),pow(MsqU[de]/Mch[ie],2.));	
	NQ11f+=Dp*temp;
	NQ21f+=Dm*temp;
	}
		NQ11f*=-4./3.*param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;

NQ21f*=4./3.*param->mass_mu*param->tan_beta*param->tan_beta/param->mass_W/param->mass_W/(param->mass_H*param->mass_H-param->mass_W*param->mass_W)*aY*mass_b_muW/sw2;


	double complex CQ1four_1=NQ11f+BQ11f;

	if(cabs(CQ1four_1)*alphas_mu/4./pi>cabs(CQ0b[1])) CQ1four_1*=cabs(CQ0b[1])/cabs(CQ1four_1)*4.*pi/alphas_mu;
		
	CQ1b[1]+=CQ1four_1;

	double complex CQ2four_1=NQ21f+BQ21f;
	
	if(cabs(CQ2four_1)*alphas_mu/4./pi>cabs(CQ0b[2])) CQ2four_1*=cabs(CQ0b[2])/cabs(CQ2four_1)*4.*pi/alphas_mu;
	CQ1b[2]+=CQ2four_1;

	CQ0b[1]*=pow(eta_mu,-4./beta0);
	CQ0b[2]*=pow(eta_mu,-4./beta0);
	CQ1b[1]*=pow(eta_mu,-4./beta0)*eta_mu;
	CQ1b[2]*=pow(eta_mu,-4./beta0)*eta_mu;
	}
	
/* NMSSM */

	if((param->mass_A02!=0.)||(param->mass_H03!=0.))
	{
		double s=param->lambdaSNMSSM/param->lambdaNMSSM;
		double v=sqrt(1./sqrt(2.)/param->Gfermi);
		
		double v_deltam_s=v/s*(sqrt(2.)*param->AlambdaNMSSM-2.*param->kappaNMSSM*s)/(sqrt(2.)*param->AlambdaNMSSM+param->kappaNMSSM*s);
		
		CQ0b[1]=0.;
		CQ0b[2]=0.;
		CQ1b[1]=0.;
		CQ1b[2]=0.;
		
		double mH0[4],mA0[3],mstop[3];
		
		mstop[0]=param->mass_upr;
		mstop[1]=param->mass_t1;
		mstop[2]=param->mass_t2;
		
		mH0[1]=param->mass_h0;
		mH0[2]=param->mass_H0;
		mH0[3]=param->mass_H03;
		mA0[1]=param->mass_A0;
		mA0[2]=param->mass_A02;
		
		double Ralj[3][3][3],Qalj[4][3][3],G1[4][4][3][3];
		double T1[3][4][4],T2[4][4][4];
		
		double TU[4][4];
		for(ie=1;ie<=3;ie++) TU[ie][1]=TU[1][ie]=0.;
		TU[1][1]=1.;
		for(ie=1;ie<=2;ie++) for(je=1;je<=2;je++) TU[ie+1][je+1]=param->stop_mix[ie][je];

		
		double vu=sqrt(pow(sin(atan(param->tan_beta)),2.)/sqrt(2.)/param->Gfermi);
		double vd=vu/param->tan_beta;
		int le;

		for(ae=1;ae<=2;ae++) for(le=1;le<=2;le++) for(je=1;je<=2;je++) Ralj[ae][le][je]=-param->g2/sqrt(2.)*(param->A0_mix[ae][1]*param->charg_Umix[2][le]*param->charg_Vmix[2][je]+param->A0_mix[ae][2]*param->charg_Umix[1][le]*param->charg_Vmix[2][je])-param->lambdaNMSSM/sqrt(2.)*param->A0_mix[ae][3]*param->charg_Umix[2][le]*param->charg_Vmix[2][je];
		
		for(ae=1;ae<=3;ae++) for(le=1;le<=2;le++) for(je=1;je<=2;je++) Qalj[ae][le][je]=param->g2/sqrt(2.)*(param->H0_mix[ae][1]*param->charg_Umix[2][le]*param->charg_Vmix[2][je]+param->H0_mix[ae][2]*param->charg_Umix[1][le]*param->charg_Vmix[2][je])-param->lambdaNMSSM/sqrt(2.)*param->H0_mix[ae][3]*param->charg_Umix[2][le]*param->charg_Vmix[2][je];
		
		for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) for(je=1;je<=2;je++) for(le=1;le<=2;le++) G1[ie][ke][je][le]=(TU[ie][2]*TU[ke][2]-kron(ie,1)*kron(ke,1))*param->charg_Vmix[1][le]*param->charg_Umix[2][je]-mass_top_muW/sqrt(2.)/sin(atan(param->tan_beta))/param->mass_W*TU[ie][3]*TU[ke][2]*param->charg_Vmix[2][le]*param->charg_Umix[2][je];
		
		for(ae=1;ae<=2;ae++) for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) T1[ae][ie][ke]=(TU[ie][3]*TU[ke][2]-TU[ie][2]*TU[ke][3])*((param->lambdaNMSSM/sqrt(2.)*(vd*param->A0_mix[ae][3]+s*param->A0_mix[ae][1]))-param->A_u*param->A0_mix[ae][2]);
		
		for(ae=1;ae<=3;ae++) for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) T2[ae][ie][ke]=-mass_top_muW/2./param->mass_W*(2.*mass_top_muW*param->H0_mix[ae][2]*(TU[ie][2]*TU[ke][2]+TU[ie][3]*TU[ke][3])	+((param->lambdaNMSSM/sqrt(2.)*(vd*param->H0_mix[ae][3]+s*param->H0_mix[ae][1]))+param->A_u*param->H0_mix[ae][2])*(TU[ie][3]*TU[ke][2]+TU[ie][2]*TU[ke][3]))	+param->mass_Z/2./sqrt(1.-sw2)*(1.-4./3.*sw2)*param->H0_mix[ae][2]*(TU[ie][1]*TU[ke][1]+TU[ie][2]*TU[ke][2])
	+2./3.*param->mass_W*sw2/(1.-sw2)*param->H0_mix[ae][2]*TU[ie][3]*TU[ke][3];
			
		double complex CQ1H=0.;
		for(ae=1;ae<=3;ae++) CQ1H+=(param->mass_H*param->mass_H/param->mass_W/param->mass_W*param->H0_mix[ae][1]*param->H0_mix[ae][1]*f30(param->mass_H*param->mass_H/mass_top_muW/mass_top_muW,param->mass_W*param->mass_W/mass_top_muW/mass_top_muW)	+mass_top_muW*mass_top_muW*mH0[ae]*mH0[ae]/param->mass_W/param->mass_W/param->mass_H/param->mass_H*f30(mass_top_muW*mass_top_muW/param->mass_H/param->mass_H,mass_top_muW*mass_top_muW/param->mass_W/param->mass_W))
		/mH0[ae]/mH0[ae];
		
		CQ1H*=-param->mass_mu/4.*param->tan_beta*param->tan_beta;
		
		double complex CQ2H=0.;
		for(ae=1;ae<=2;ae++) CQ2H+=((param->mass_H*param->mass_H/param->mass_W/param->mass_W*param->A0_mix[ae][1]*param->A0_mix[ae][1]+kron(ae,2)*param->A0_mix[ae][1])*f30(param->mass_H*param->mass_H/mass_top_muW/mass_top_muW,param->mass_W*param->mass_W/mass_top_muW/mass_top_muW)	+mass_top_muW*mass_top_muW*mA0[ae]*mA0[ae]/param->mass_W/param->mass_W/param->mass_H/param->mass_H*f30(mass_top_muW*mass_top_muW/param->mass_H/param->mass_H,mass_top_muW*mass_top_muW/param->mass_W/param->mass_W))/mA0[ae]/mA0[ae];
		
		CQ2H*=param->mass_mu/4.*param->tan_beta*param->tan_beta;
		
		double complex CAH=-I*param->lambdaNMSSM*param->AlambdaNMSSM/param->g2/param->mass_W*param->tan_beta*f30(param->mass_H*param->mass_H/mass_top_muW/mass_top_muW,param->mass_W*param->mass_W/mass_top_muW/mass_top_muW);
		
		double complex CQ1c=0.;
		
		for(ae=1;ae<=3;ae++) for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) for(je=1;je<=2;je++) for(le=1;le<=2;le++) CQ1c+=G1[ie][ke][je][le]/mH0[ae]/mH0[ae]*( sqrt(2.)*param->H0_mix[ae][1]*param->H0_mix[ae][1]*Mch[je]/param->mass_W/cos(atan(param->tan_beta))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/Mch[je],2.))
			-2.*sqrt(2.)*param->H0_mix[ae][1]/param->g2*kron(ie,ke)*(Qalj[ae][le][je]*f40(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.))+Mch[je]/Mch[le]*Qalj[ae][je][le]*f30(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.)))		+2.*sqrt(2.)*param->H0_mix[ae][1]*T2[ae][ie][ke]*Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow(Mch[je]/mstop[ke-1],2.))
			+mH0[ae]*mH0[ae]/Mch[je]/Mch[je]*kron(ie,ke)*(param->charg_Umix[2][je]*param->charg_Vmix[1][le]*f50(pow(mstop[ie-1]/Mch[je],2.),pow(Mch[le]/Mch[je],2.),pow(param->mass_nutl/Mch[le],2.))
		
		-Mch[le]/Mch[je]*param->charg_Umix[2][le]*param->charg_Vmix[1][je]* f60(pow(mstop[ie-1]/Mch[je],2.),pow(Mch[le]/Mch[je],2.),pow(param->mass_nutl/Mch[le],2.))));
			
		CQ1c*=param->mass_mu/4.*param->tan_beta*param->tan_beta;
		
		double complex CQ2c=0.;
		
		for(ae=1;ae<=2;ae++) for(ie=1;ie<=3;ie++) for(ke=1;ke<=3;ke++) for(je=1;je<=2;je++) for(le=1;le<=2;le++) CQ2c+=G1[ie][ke][je][le]/mA0[ae]/mA0[ae]*(sqrt(2.)*param->A0_mix[ae][1]*param->A0_mix[ae][1]*Mch[je]/param->mass_W/cos(atan(param->tan_beta))*kron(ie,ke)*kron(le,je)*f80(pow(mstop[ie-1]/Mch[je],2.))
			-2.*sqrt(2.)*param->A0_mix[ae][1]/param->g2*kron(ie,ke)*(-Ralj[ae][le][je]*f40(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.))+Mch[je]/Mch[le]*Ralj[ae][je][le]*f30(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.)))			-sqrt(2.)*param->A0_mix[ae][1]*T1[ae][ie][ke]*mass_top_muW*Mch[je]/mstop[ke-1]/mstop[ke-1]*kron(le,je)*f30(pow(mstop[ie-1]/mstop[ke-1],2.),pow(Mch[je]/mstop[ke-1],2.))
				+mA0[ae]*mA0[ae]/Mch[je]/Mch[je]*kron(ie,ke)*(param->charg_Umix[2][je]*param->charg_Vmix[1][le]*f50(pow(mstop[ie-1]/Mch[je],2.),pow(Mch[le]/Mch[je],2.),pow(param->mass_nutl/Mch[le],2.))
			-Mch[le]/Mch[je]*param->charg_Umix[2][le]*param->charg_Vmix[1][je]*f60(pow(mstop[ie-1]/Mch[je],2.),pow(Mch[le]/Mch[je],2.),pow(param->mass_nutl/Mch[le],2.))));
			
		CQ2c*=-param->mass_mu/4.*param->tan_beta*param->tan_beta;		
				
		double complex CAc=0.;
		for(ie=1;ie<=3;ie++) for(je=1;je<=2;je++) for(le=1;le<=2;le++) CAc+=I*param->tan_beta/sqrt(2.)*G1[ie][ie][je][le]*(v_deltam_s*kron(le,je)*fabs(Mch[je]/param->mass_W)*f80(pow(mstop[ie-1]/Mch[je],2.))-(Ralj[1][je][le]*fabs(Mch[je]/Mch[le])*f30(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.))-Ralj[1][le][je]*f40(pow(mstop[ie-1]/Mch[le],2.),pow(Mch[je]/Mch[le],2.))));
				
		CQ0b[1]=(CQ1H+CQ1c)*mass_b_muW/sw2/epsfac;
		CQ0b[2]=(CQ2H+CQ2c)*mass_b_muW/sw2/epsfac;

		double complex CA=CAH+CAc;

		if(param->mass_A0>mu_W) CQ0b[2]+=-v_deltam_s/2.*mass_b_muW/sw2*param->mass_mu*CA/param->mass_A0/param->mass_A0;

		CQ0b[1]*=pow(eta_mu,-4./beta0);
		CQ0b[2]*=pow(eta_mu,-4./beta0);
		CQ1b[1]*=pow(eta_mu,-4./beta0)*eta_mu;
		CQ1b[2]*=pow(eta_mu,-4./beta0)*eta_mu;
		
		if((param->mass_A0>param->mass_b_pole)&&(param->mass_A0<mu_W))
		{	
			double alphas_Ma1=alphas_running(param->mass_A0,param->mass_top_pole,param->mass_b_pole,param);	
			double eta_a1=alphas_Ma1/alphas_mu;
			double mass_b_ma1=running_mass(param->mass_b,param->mass_b,param->mass_A0,param->mass_top_pole,param->mass_b,param);
			CQ0b[2]+=-v_deltam_s/2.*mass_b_ma1/sw2*param->mass_mu*CA/param->mass_A0/param->mass_A0*pow(eta_a1,-4./beta0);
		}
		
		if(param->mass_A0<param->mass_b_pole)
		{	
			double width_A0=1.e-6;	
CQ0b[2]+=v_deltam_s/2.*param->mass_b/sw2*param->mass_mu*CA/(param->m_Bs*param->m_Bs-param->mass_A0*param->mass_A0+I*param->mass_A0*width_A0);
		}
	}
			
	return;
}