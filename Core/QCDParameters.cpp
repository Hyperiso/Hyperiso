#include "QCDParameters.h"

QCDParameters::QCDParameters(double alpha_Z, double m_Z, double masst_pole, double massb_b, double mass_u, double mass_d, double mass_s, double mass_c) {
    
    this->mass_Z = m_Z;
    this->alphas_MZ = alpha_Z;
    this->mass_t_pole = masst_pole;
    this->mass_b_b = massb_b;
    this->mass_c = mass_c;
    this->mass_s = mass_s;

    Logger* logger = Logger::getInstance();

    logger->info("m_Z : " + std::to_string(m_Z));
    logger->info("alpha_MZ : " + std::to_string(alpha_Z));
    logger->info("mass_c : " + std::to_string(mass_c));
    logger->info("mass_s : " + std::to_string(mass_s));


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

    Logger* logger = Logger::getInstance();

    double mass_b, mass_t;

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
        if (Lambda5 == -1){
        logger->error("Error for Lambda5 calculation in QCCDParameters.cpp");
            return Lambda5;}
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

    logger->debug("Invalid parameters or conditions in alphas_running function");
}

std::tuple<double, double, double> calculateBetas(int nf) {
    double beta0 = 11.0 - 2.0 / 3.0 * nf;
    double beta1 = 51.0 - 19.0 / 3.0 * nf;
    double beta2 = 2857.0 - 5033.0 / 9.0 * nf + 325.0 / 27.0 * nf * nf;
    return {beta0, beta1, beta2};
}

std::tuple<double, double, double> calculateGammas(int nf) {
    double gamma0 = 2.0;
    double gamma1 = 101.0 / 12.0 - 5.0 / 18.0 * nf;
    double gamma2 = 1.0 / 32.0 * (1249.0 - (2216.0 / 27.0 + 160.0 / 3.0 * ZETA3) * nf - 140.0 / 81.0 * nf * nf);
    return {gamma0, gamma1, gamma2};
}

double calculateR(double alphas, double beta0, double beta1, double beta2, double gamma0, double gamma1, double gamma2) {
    double term1 = pow(beta0 / (2.0 * PI) * alphas, 2.0 * gamma0 / beta0);
    double term2 = 1.0 + (2.0 * gamma1 / beta0 - beta1 * gamma0 / (beta0 * beta0)) * alphas / PI;
    double term3 = 0.5 * (pow(2.0 * gamma1 / beta0 - beta1 * gamma0 / (beta0 * beta0), 2.0)
                        + 2.0 * gamma2 / beta0
                        - beta1 * gamma1 / (beta0 * beta0)
                        - beta2 * gamma0 / (16.0 * beta0 * beta0)
                        + beta1 * beta1 * gamma0 / (2.0 * beta0 * beta0 * beta0)) * pow(alphas / PI, 2.0);
    return term1 * (term2 + term3);
}

double QCDParameters::running_mass(double quark_mass, double Qinit, double Qfin,  std::string option_massb, std::string option_masst)
/* computes the running quark mass at the energy Qfin, from a given running quark mass quark_mass at energy Qinit */
{

    double mbot = (option_massb == "pole") ? mass_b_pole : mass_b_b;
    double mtop = (option_masst == "pole") ? mass_t_pole : mass_t_t;

    double alphas_Qinit = runningAlphasCalculation(Qinit, option_massb, option_masst);

    int nf_init = (Qinit <= mbot) ? 4 : (Qinit <= mtop) ? 5 : 6;
    int nf_fin = (Qfin <= mbot) ? 4 : (Qfin <= mtop) ? 5 : 6;

    auto [beta0_init, beta1_init, beta2_init] = calculateBetas(nf_init);
    auto [gamma0_init, gamma1_init, gamma2_init] = calculateGammas(nf_init);

    double R_Qinit = calculateR(alphas_Qinit, beta0_init, beta1_init, beta2_init, gamma0_init, gamma1_init, gamma2_init);

    double intermediate_mass = quark_mass;

    auto runStep = [&](double Qstart, double Qend, int nf_start, int nf_end, double &mass) {
        double alphas_start = runningAlphasCalculation(Qstart, option_massb, option_masst);
        double alphas_end = runningAlphasCalculation(Qend, option_massb, option_masst);

        auto [beta0_start, beta1_start, beta2_start] = calculateBetas(nf_start);
        auto [gamma0_start, gamma1_start, gamma2_start] = calculateGammas(nf_start);
        double R_start = calculateR(alphas_start, beta0_start, beta1_start, beta2_start, gamma0_start, gamma1_start, gamma2_start);

        auto [beta0_end, beta1_end, beta2_end] = calculateBetas(nf_end);
        auto [gamma0_end, gamma1_end, gamma2_end] = calculateGammas(nf_end);
        double R_end = calculateR(alphas_end, beta0_end, beta1_end, beta2_end, gamma0_end, gamma1_end, gamma2_end);

        mass *= R_end / R_start;
    };

    if (nf_init != nf_fin) {
        if (nf_init == 4 && nf_fin == 6) {
            runStep(Qinit, mbot, 4, 5, intermediate_mass);
            runStep(mbot, mtop, 5, 6, intermediate_mass);
        } else if (nf_init == 6 && nf_fin == 4) {
            runStep(Qinit, mtop, 6, 5, intermediate_mass);
            runStep(mtop, mbot, 5, 4, intermediate_mass);
        } else if (nf_init < nf_fin) {
            double Q_intermediate = (nf_init == 4) ? mbot : mtop;
            runStep(Qinit, Q_intermediate, nf_init, nf_init + 1, intermediate_mass);
        } else {
            double Q_intermediate = (nf_fin == 4) ? mbot : mtop;
            runStep(Qinit, Q_intermediate, nf_init, nf_init - 1, intermediate_mass);
        }
    }

    if (nf_init != nf_fin) {
        runStep((nf_init < nf_fin) ? ((nf_init == 4) ? mbot : mtop) : ((nf_fin == 4) ? mbot : mtop), Qfin, nf_fin, nf_fin, intermediate_mass);
    } else {
        runStep(Qinit, Qfin, nf_init, nf_fin, intermediate_mass);
    }

    return intermediate_mass;

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
	double a2=(4343./162.+4.*PI*PI-pow(PI,4.)/4.+22./3.*ZETA3)*9.-(1798./81.+56./3.*ZETA3)*6.
	-(55./3.-16.*ZETA3)*8./3.+1600./81.;
	
	return mb_pole*(1.-2./9.*pow(as,2.)
	/* -2./9.*pow(as,3.)/PI*(beta0*(L+1.)+a1/2.) */
	)
	-2./9.*pow(as,2.)
	/* -2./9.*pow(as,4.)/PI/PI*(beta0*beta0*(3./4.*L*L+L+ZETA3/2.+PI*PI/24.+1./4.)
	+beta0*a1/2.*(3./2.*L+1.)+beta1/4.*(L+1.)+a1*a1/16.+a2/8.+(3.-1./36.)*4./3.*PI*PI) */
	;
}

/*--------------------------------------------------------------------*/

double QCDParameters::mt_mt(double mt_pole)
/* computes the top mass mt(mt) */
{
	double alphas_mtop=runningAlphasCalculation(mt_pole, "running"); 
    this->mass_t_t = this->mass_t_pole/(1.+alphas_mtop/PI*(4./3.+alphas_mtop/PI*(307./32.+PI*PI/3.+PI*PI/9.*log(2.)-1./6.*ZETA3-71./144.*5.)));
    alphas_mtop=runningAlphasCalculation(this->mass_t_t,"running","running");
    this->mass_t_t = this->mass_t_pole/(1.+alphas_mtop/PI*(4./3.+alphas_mtop/PI*(307./32.+PI*PI/3.+PI*PI/9.*log(2.)-1./6.*ZETA3-71./144.*5.)));
	return mass_t_t;

}

