#include "./QCDParameters.h"
#include <iostream>

QCDParameters::QCDParameters(double alpha_Z, double m_Z, double masst_pole, double massb_b, double mass_u, double mass_d, double mass_s, double mass_c) {
    mass_Z = m_Z;
    alphas_MZ = alpha_Z;
    // mass_b_pole=massb_pole;
    mass_t_pole = masst_pole;
    mass_b_b = massb_b;
    // mass_t_t = masst_t;
    this->mt_mt(this->mass_t_pole);
    this->mb_pole(mass_b_b, mass_u, mass_d, mass_s, mass_c);
}

double QCDParameters::alphasRunning(double Q, double Lambda, int nf) const{
    double beta0 = 11.-2./3.*nf;
    double beta1=51.-19./3.*nf;
    double beta2=2857.-5033.*nf/9.+325./27.*nf*nf;

    return 4.*PI/beta0/log(pow(Q/Lambda,2.))*(1.-2.*beta1/beta0/beta0*log(log(pow(Q/Lambda,2.)))/log(pow(Q/Lambda,2.))+4.*beta1*beta1/pow(beta0*beta0*log(pow(Q/Lambda,2.)),2.)*(pow(log(log(pow(Q/Lambda,2.)))-1./2.,2.)+beta2*beta0/8./beta1/beta1-5./4.));
}

double QCDParameters::matchLambda(double alpha_running, double Q, int nf){

    double Lambda_min=1.e-3;
    double Lambda_max=1.;
    double alphas_min=0.;
    double alphas_max=0;
    double alphas_moy = 0;
    double Lambda_moy = 0;
    while((fabs(1.-alphas_min/alpha_running)>=1.e-4)&&(fabs(1.-Lambda_min/Lambda_max)>1.e-5)) {
        alphas_min = alphasRunning(Q,Lambda_min, nf);
        alphas_max = alphasRunning(Q,Lambda_max, nf);

        Lambda_moy=(Lambda_min+Lambda_max)/2.;
        alphas_moy = alphasRunning(Q,Lambda_moy, nf);

        if((alpha_running>=alphas_min)&&(alpha_running<=alphas_moy))
            Lambda_max=Lambda_moy;
        else Lambda_min=Lambda_moy;
    }
    
    if(fabs(1.-Lambda_min/Lambda_max)<=1.e-5) {
        Lambda5=-1.;
        return -1.;
    }
    return Lambda_min;
}

double QCDParameters::runningAlphasCalculation(double Q, std::string option_massb, std::string option_masst){

    if (option_massb == "pole")
        mass_b = mass_b_pole;
    else if (option_massb == "running")
        mass_b = mass_b_b;

    if (option_masst == "pole")
        mass_t = mass_t_pole;
    else if (option_masst == "running")
        mass_t = mass_t_t;

    if (Lambda5 == -1.0) {
        return -1.0;
    }

    if (Lambda5 == 0.0) {
        Lambda5 = matchLambda(alphas_MZ, mass_Z, 5);
        if (Lambda5 == -1);
            return Lambda5;
    }
    double alphas_running = alphasRunning(Q, Lambda5, 5);

    if (Q <= mass_t && Q >= mass_b) {
        // 5 active flavors
        return alphas_running;
    } else if (Q > mass_t) {
        // 6 active flavors
        nf = 6;
        alphas_running = alphasRunning(mass_t, Lambda5, 5);
        Lambda6 = matchLambda( alphas_running,  mass_t, nf);
        return alphasRunning(Q,Lambda6, nf);
    } else {
        // 4 active flavors
        nf = 4;
        alphas_running = alphasRunning(mass_b, Lambda5, 5);
        Lambda4 = matchLambda(alphas_running, mass_b, nf);
        return alphasRunning(Q, Lambda4, nf);
    }

    throw std::logic_error("Invalid parameters or conditions in alphas_running function");
}


#define zeta3 1.2020569031595942855
#define PI    3.1415926535897932385

double QCDParameters::running_mass(double quark_mass, double Qinit, double Qfin,  std::string option_massb, std::string option_masst)
/* computes the running quark mass at the energy Qfin, from a given running quark mass quark_mass at energy Qinit */
{
    double mbot, mtop;

    if (option_massb == "pole")
        mbot = mass_b_pole;
    else if (option_massb == "running")
        mbot = mass_b_b;

    if (option_masst == "pole")
        mtop = mass_t_pole;
    else if (option_masst == "running")
        mtop = mass_t_t;
	double alphas_Qinit,alphas_Qfin,running_mass;
	double beta0_init,beta1_init,beta2_init,gamma0_init,gamma1_init,gamma2_init;
    double beta0_fin,beta1_fin,beta2_fin,gamma0_fin,gamma1_fin,gamma2_fin;
	int nf_init, nf_fin;

	alphas_Qinit=runningAlphasCalculation(Qinit, option_massb, option_masst);
	
	double R_Qinit,R_Qfin;
	
    if (Qinit <= mbot) {
        nf_init=4;
    } else if(Qinit <= mtop) {
        nf_init = 5;
    } else {
        nf_init = 6;
    }

    if (Qfin <= mbot) {
        nf_fin=4;
    } else if(Qfin <= mtop) {
        nf_fin = 5;
    } else {
        nf_fin = 6;
    }

    beta0_init=11.-2./3.*nf_init;
    beta1_init=51.-19./3.*nf_init;
    beta2_init=2857.-5033./9.*nf_init+325./27.*nf_init*nf_init;

    gamma0_init=2.;
    gamma1_init=101./12.-5./18.*nf_init;
    gamma2_init=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf_init-140./81.*nf_init*nf_init);

    beta0_fin=11.-2./3.*nf_fin;
    beta1_fin=51.-19./3.*nf_fin;
    beta2_fin=2857.-5033./9.*nf_fin+325./27.*nf_fin*nf_fin;

    gamma0_fin=2.;
    gamma1_fin=101./12.-5./18.*nf_fin;
    gamma2_fin=1./32.*(1249.-(2216./27.+160./3.*zeta3)*nf_fin-140./81.*nf_fin*nf_fin);

    R_Qinit = pow(beta0_init/2./PI*alphas_Qinit,2.*gamma0_init/beta0_init)*(1.+(2.*gamma1_init/beta0_init-beta1_init*gamma0_init/(beta0_init*beta0_init))*alphas_Qinit/PI+1./2.*(pow(2.*gamma1_init/beta0_init-beta1_init*gamma0_init/(beta0_init*beta0_init),2.)+2.*gamma2_init/beta0_init-beta1_init*gamma1_init/(beta0_init*beta0_init)-beta2_init*gamma0_init/16./beta0_init/beta0_init+beta1_init*beta1_init*gamma0_init/2./(beta0_init*beta0_init*beta0_init))*pow(alphas_Qinit/PI,2.));
	
    alphas_Qfin=runningAlphasCalculation(Qfin, option_massb, option_masst);

    R_Qfin = pow(beta0_fin/2./PI*alphas_Qfin,2.*gamma0_fin/beta0_fin)*(1.+(2.*gamma1_fin/beta0_fin-beta1_fin*gamma0_fin/(beta0_fin*beta0_fin))*alphas_Qfin/PI+1./2.*(pow(2.*gamma1_fin/beta0_fin-beta1_fin*gamma0_fin/(beta0_fin*beta0_fin),2.)+2.*gamma2_init/beta0_init-beta1_fin*gamma1_fin/(beta0_fin*beta0_fin)-beta2_fin*gamma0_fin/16./beta0_fin/beta0_fin+beta1_fin*beta1_fin*gamma0_fin/2./(beta0_fin*beta0_fin*beta0_fin))*pow(alphas_Qfin/PI,2.));
	
    return R_Qfin/R_Qinit*quark_mass;

}

/*--------------------------------------------------------------------*/

double QCDParameters::mb_pole(double mass_b, double mass_u, double mass_d, double mass_s, double mass_c)
/* computes the b pole mass */
{
	double alphas_mb=runningAlphasCalculation(mass_b, "running", "pole");	
    this->mass_b_pole = mass_b*(1.+alphas_mb/PI*(4./3.+alphas_mb/PI*((13.4434-1.0414*4.+1.0414*4./3.*((mass_u+mass_d+mass_s+mass_c)/mass_b)))));
 	return this->mass_b_pole;
	/* +alphas_mb/PI*(190.595-4.*(26.655-4.*0.6527)) */
}

/*--------------------------------------------------------------------*/

double QCDParameters::mc_pole(double mass_u, double mass_d, double mass_s, double mass_c)
/* computes the c pole mass */{
	double alphas_mc=runningAlphasCalculation(mass_c);	

 	return mass_c*(1+alphas_mc/PI*(4./3.+alphas_mc/PI*((13.4434-1.0414*3.+1.0414*4./3.*((mass_u+mass_d+mass_s)/mass_c)))));
}

double QCDParameters::mb_1S(double mb_pole)
/* computes the 1S b mass */
{
	double mu=mb_pole/2.;
	double as=runningAlphasCalculation(mu);
	double L=log(mu/(4./3.*as*mb_pole));
	double beta0=11.-8./3.;
	double beta1=62.-32./3.;
	double a1=31./3.-40./9.;
	double a2=(4343./162.+4.*PI*PI-pow(PI,4.)/4.+22./3.*zeta3)*9.-(1798./81.+56./3.*zeta3)*6.
	-(55./3.-16.*zeta3)*8./3.+1600./81.;
	
	return mb_pole*(1.-2./9.*pow(as,2.)
	/* -2./9.*pow(as,3.)/PI*(beta0*(L+1.)+a1/2.) */
	)
	-2./9.*pow(as,2.)
	/* -2./9.*pow(as,4.)/PI/PI*(beta0*beta0*(3./4.*L*L+L+zeta3/2.+PI*PI/24.+1./4.)
	+beta0*a1/2.*(3./2.*L+1.)+beta1/4.*(L+1.)+a1*a1/16.+a2/8.+(3.-1./36.)*4./3.*PI*PI) */
	;
}

/*--------------------------------------------------------------------*/

double QCDParameters::mt_mt(double mt_pole)
/* computes the top mass mt(mt) */
{
	double alphas_mtop=runningAlphasCalculation(mt_pole, "running"); 
    this->mass_t_t = this->mass_t_pole/(1.+alphas_mtop/PI*(4./3.+alphas_mtop/PI*(307./32.+PI*PI/3.+PI*PI/9.*log(2.)-1./6.*zeta3-71./144.*5.)));
    alphas_mtop=runningAlphasCalculation(this->mass_t_t,"running","running");
    this->mass_t_t = this->mass_t_pole/(1.+alphas_mtop/PI*(4./3.+alphas_mtop/PI*(307./32.+PI*PI/3.+PI*PI/9.*log(2.)-1./6.*zeta3-71./144.*5.)));
	return mass_t_t;

}
