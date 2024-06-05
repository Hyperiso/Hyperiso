#include "QCDParameters.h"

QCDParameters::QCDParameters(double alpha_Z, double m_Z, double masst_pole, double massb_b, double mass_u, double mass_d, double mass_s, double mass_c) {
    
    this->Lambda5 = this->matchLambda(alpha_Z, m_Z, 5);
    this->mass_u = mass_u;
    this->mass_d = mass_d;
    this->mass_s = mass_s;
    this->mass_c = mass_c;
    this->mass_b_b = massb_b;
    this->mass_t_pole = masst_pole;
    this->mt_mt(this->mass_t_pole);
    this->mb_pole(mass_b_b, mass_u, mass_d, mass_s, mass_c);
    this->setMassTypes("pole", "pole");
}

double QCDParameters::runningAlphasCalculation(double Q, std::string option_massb, std::string option_masst) {
    Logger* logger = Logger::getInstance();
    this->setMassTypes(option_massb, option_masst);

    int n_i = 5;
    int n_f = this->getNf(Q);

    if (n_f < 4) {
        logger->warn("Scale for alpha_s calculation is below charm mass.");
    }

    double L = this->Lambda5;
    auto Q_bounds = this->getOrderedMasses();

    while (n_i > n_f) {
        double alpha_match = this->alphasRunning(Q_bounds.at(n_i - 1), L, n_i);
        L = this->matchLambda(alpha_match, Q_bounds.at(n_i - 1), n_i - 1);
        --n_i;
    }

    while (n_i < n_f) {
        double alpha_match = this->alphasRunning(Q_bounds.at(n_i), L, n_i);
        L = this->matchLambda(alpha_match, Q_bounds.at(n_i), n_i + 1);
        ++n_i;
    }

    return this->alphasRunning(Q, L, n_f);
}

double QCDParameters::running_mass(double quark_mass, double Qinit, double Qfin, std::string option_massb, std::string option_masst)
/* computes the running quark mass at the energy Qfin, from a given running quark mass quark_mass at energy Qinit */
{
    this->setMassTypes(option_massb, option_masst);
    int n_i = this->getNf(Qinit);
    int n_f = this->getNf(Qfin);
    auto Q_bounds = this->getOrderedMasses();

    while (n_i > n_f) {
        quark_mass = this->runMass(quark_mass, Qinit, Q_bounds.at(n_i - 1), n_i);
        Qinit = Q_bounds.at(n_i - 1);
        --n_i;
    }

    while (n_i < n_f) {
        quark_mass = this->runMass(quark_mass, Qinit, Q_bounds.at(n_i), n_i);
        Qinit = Q_bounds.at(n_i);
        ++n_i;
    }

    return this->runMass(quark_mass, Qinit, Qfin, n_f);
}

void QCDParameters::setMassTypes(std::string m_b_type, std::string m_t_type) {
    if (m_b_type != "")
        this->m_b_type = m_b_type;
    if (m_t_type != "")
        this->m_t_type = m_t_type;
}

std::vector<double> QCDParameters::getOrderedMasses() {
    double m_b = this->m_b_type == "running" ? this->get_mb_mb() : this-> get_mb_pole();
    double m_t = this->m_t_type == "running" ? this->get_mt_mt() : this-> get_mt_pole();
    return {this->mass_u, this->mass_d, this->mass_s, this->mass_c, m_b, m_t};
}

double QCDParameters::runMass(double mass, double Q_i, double Q_f, int nf) {
    return mass * this->R(this->runningAlphasCalculation(Q_f), nf) / this->R(this->runningAlphasCalculation(Q_i), nf);
}

int QCDParameters::getNf(double Q) {
    auto masses = this->getOrderedMasses();
    for (int i = 0; i < masses.size(); ++i) {
        if (Q < masses.at(i))
            return i;
    }
    return 6;
}

std::tuple<double, double, double> QCDParameters::getBetas(int nf) const {
    double b0 = 11 - 2. * nf / 3;
    double b1 = 51 - 19. * nf / 3;
    double b2 = 2857. - 5033. * nf / 9 + 325. * nf * nf / 27;
    return {b0, b1, b2};
}

std::tuple<double, double, double> QCDParameters::getGammas(int nf) const {
    double g0 = 2;
    double g1 = 101.0 / 12 - 5. * nf / 18;
    double g2 = (1249. - (2216. / 27 + 160 * ZETA3 / 3) * nf - 140. * nf * nf / 81) / 32;
    return {g0, g1, g2};
}

double QCDParameters::R(double alpha, int nf) const {
    auto [b0, b1, b2] = this->getBetas(nf);
    auto [g0, g1, g2] = this->getGammas(nf);
    double b02 = b0 * b0;
    double a = std::pow(b0 * alpha / (2 * PI), 2 * g0 / b0);
    double b = (2 * g1 / b0 - b1 * g0 / b02) * alpha / PI;
    double c = .5 * (std::pow(2 * g1 / b0 - b1 * g0 / b02, 2) 
                   + 2 * g2 / b0 - b1 * g1 / b02 - b2 * g0 / (16 * b02) 
                   + b1 * b1 * g0 / (2 * b0 * b02)) * std::pow(alpha / PI, 2);
    return a * (1 + b + c);
}

double QCDParameters::alphasRunning(double Q, double Lambda, int nf) const {
    double r = std::pow(Q / Lambda, 2);
    double L = std::log(r);
    double LL = std::log(L);
    auto [b0, b1, b2] = this->getBetas(nf);
    double b02 = b0 * b0;
    double b12 = b1 * b1;
    return 4 * PI * (1 - 2 * b1 * LL / (b02 * L) + 4 * b12 * (std::pow(LL - .5, 2) + b2 * b0 / 8 / b12 - 1.25) / std::pow(b02 * L, 2)) / (b0 * L);
}

double QCDParameters::matchLambda(double target_alpha, double Q, int nf){
    auto f = [&](double L) { return this->alphasRunning(Q, L, nf) - target_alpha; };
    double L_min = 1e-3;
    double L_max = 1;
    double L_moy = L_min;

    while (std::abs(1 - L_min / L_max) > 1e-5) {
        L_moy = (L_min + L_max) / 2;
        f(L_moy) > 0 ? L_max = L_moy : L_min = L_moy;
    }

    if (std::abs(f(L_moy)) > 1e-5) {
        Logger::getInstance()->error("Unable to find suitable QCD Lambda value to match alpha_s = " + std::to_string(target_alpha) 
                 + " at scale " + std::to_string(Q) + " GeV with " + std::to_string(nf) + " active flavors.");
        return -1;
    }

    return L_min;
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
