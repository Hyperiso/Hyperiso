#include "QCDParameters.h"

#define zeta3 1.2020569031595942855
#define pi    3.1415926535897932385

double running_mass(double quark_mass, double Qinit, double Qfin,  double mtop, double mbot, struct parameters* param)
/* computes the running quark mass at the energy Qfin, from a given running quark mass quark_mass at energy Qinit */
{
    QCDParameters run = QCDParameters();
	double alphas_Qinit,alphas_Qfin,running_mass;
	double beta0_init,beta1_init,beta2_init,gamma0_init,gamma1_init,gamma2_init;
    double beta0_fin,beta1_fin,beta2_fin,gamma0_fin,gamma1_fin,gamma2_fin;
	int nf_init, nf_fin;

	alphas_Qinit=run.runningAlphasCalculation(Qinit);
	
	double R_Qinit,R_Qfin;
	
    if (Qinit <= mbot) {
        nf_init=4;
    }
    else if(Qinit <= mtop) {
        nf_init = 5;
    }
    else {
        nf_init = 6;
    }
    if (Qinit <= mbot) {
        nf_fin=4;
    }
    else if(Qinit <= mtop) {
        nf_fin = 5;
    }
    else {
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

    R_Qinit = pow(beta0_init/2./pi*alphas_Qinit,2.*gamma0_init/beta0_init)*(1.+(2.*gamma1_init/beta0_init-beta1_init*gamma0_init/(beta0_init*beta0_init))*alphas_Qinit/pi+1./2.*(pow(2.*gamma1_init/beta0_init-beta1_init*gamma0_init/(beta0_init*beta0_init),2.)+2.*gamma2_init/beta0_init-beta1_init*gamma1_init/(beta0_init*beta0_init)-beta2_init*gamma0_init/16./beta0_init/beta0_init+beta1_init*beta1_init*gamma0_init/2./(beta0_init*beta0_init*beta0_init))*pow(alphas_Qinit/pi,2.));
	
    alphas_Qfin=run.runningAlphasCalculation(Qfin);

    R_Qfin = pow(beta0_fin/2./pi*alphas_Qfin,2.*gamma0_fin/beta0_fin)*(1.+(2.*gamma1_fin/beta0_fin-beta1_fin*gamma0_fin/(beta0_fin*beta0_fin))*alphas_Qfin/pi+1./2.*(pow(2.*gamma1_fin/beta0_fin-beta1_fin*gamma0_fin/(beta0_fin*beta0_fin),2.)+2.*gamma2_init/beta0_init-beta1_fin*gamma1_fin/(beta0_fin*beta0_fin)-beta2_fin*gamma0_fin/16./beta0_fin/beta0_fin+beta1_fin*beta1_fin*gamma0_fin/2./(beta0_fin*beta0_fin*beta0_fin))*pow(alphas_Qfin/pi,2.));
	
    return R_Qfin/R_Qinit*quark_mass;

}

/*--------------------------------------------------------------------*/

double mb_pole(double mass_b, double mass_u, double mass_d, double mass_s, double mass_c)
/* computes the b pole mass */
{
	QCDParameters run = QCDParameters();
	double alphas_mb=run.runningAlphasCalculation(mass_b);	

 	return mass_b*(1.+alphas_mb/pi*(4./3.+alphas_mb/pi*((13.4434-1.0414*4.+1.0414*4./3.*((mass_u+mass_d+mass_s+mass_c)/mass_b))
	/* +alphas_mb/pi*(190.595-4.*(26.655-4.*0.6527)) */
	)));
}

/*--------------------------------------------------------------------*/

double mc_pole(double mass_u, double mass_d, double mass_s, double mass_c)
/* computes the c pole mass */{
	QCDParameters run = QCDParameters();
	double alphas_mc=run.runningAlphasCalculation(mass_c);	

 	return mass_c*(1+alphas_mc/pi*(4./3.+alphas_mc/pi*((13.4434-1.0414*3.+1.0414*4./3.*((mass_u+mass_d+mass_s)/mass_c)))));
}

double mb_1S(double mb_pole)
/* computes the 1S b mass */
{
    QCDParameters run = QCDParameters();
	double mu=mb_pole/2.;
	double as=run.runningAlphasCalculation(mu);
	double L=log(mu/(4./3.*as*mb_pole));
	double beta0=11.-8./3.;
	double beta1=62.-32./3.;
	double a1=31./3.-40./9.;
	double a2=(4343./162.+4.*pi*pi-pow(pi,4.)/4.+22./3.*zeta3)*9.-(1798./81.+56./3.*zeta3)*6.
	-(55./3.-16.*zeta3)*8./3.+1600./81.;
	
	return mb_pole*(1.-2./9.*pow(as,2.)
	/* -2./9.*pow(as,3.)/pi*(beta0*(L+1.)+a1/2.) */
	)
	-2./9.*pow(as,2.)
	/* -2./9.*pow(as,4.)/pi/pi*(beta0*beta0*(3./4.*L*L+L+zeta3/2.+pi*pi/24.+1./4.)
	+beta0*a1/2.*(3./2.*L+1.)+beta1/4.*(L+1.)+a1*a1/16.+a2/8.+(3.-1./36.)*4./3.*pi*pi) */
	;
}

/*--------------------------------------------------------------------*/

double mt_mt(double mt_pole)
/* computes the top mass mt(mt) */
{
	QCDParameters run = QCDParameters();
	double alphas_mtop=run.runningAlphasCalculation(mt_pole); 

 	return mt_pole*(1.-4./3.*alphas_mtop/pi);
}