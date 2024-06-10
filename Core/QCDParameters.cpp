#include "QCDParameters.h"

/**
 * @brief Creates a fully-initialized QCD Runner.
 * 
 * @param alpha_Z Value of alpha_s at M_Z in MSbar scheme 
 * @param m_Z Mass of the Z boson
 * @param masst_pole Pole mass of the top quark
 * @param massb_b Running mass of the bottom quark at m_b in MSbar
 * @param mass_u Running mass of the up quark at 2 GeV in MSbar
 * @param mass_d Running mass of the down quark at 2 GeV in MSbar
 * @param mass_s Running mass of the strange quark at 2 GeV in MSbar
 * @param mass_c Running mass of the charm quark at m_c in MSbar
 */ 
QCDParameters::QCDParameters(double alpha_Z, double m_Z, double masst_pole, double massb_b, double mass_u, double mass_d, double mass_s, double mass_c) {
    this->Lambda5 = this->matchLambda(alpha_Z, m_Z, 5);
    this->mass_u = mass_u;
    this->mass_d = mass_d;
    this->mass_s = mass_s;
    this->mass_c = mass_c;
    this->mass_b_b = massb_b;
    this->mass_t_pole = masst_pole;
    this->mt_mt();
    this->mb_pole();
    this->setMassTypes("pole", "pole");

    Logger::getInstance()->info("In QCDParameters constructor mb(81 GeV) = " + std::to_string(this->running_mass(4.25, 4.25, 81, "running")));
}

/**
 * @brief Computes the value of alpha_s at scale Q in the MSbar scheme. 
 * 
 * @param Q Energy scale
 * @param option_massb Toggle to choose whether to take the running or pole bottom mass
 * @param option_masst Toggle to choose whether to take the running or pole top mass
 * @return The value of alpha_s(Q) in MSbar
 */
double QCDParameters::runningAlphasCalculation(double Q, std::string option_massb, std::string option_masst) {
    this->setMassTypes(option_massb, option_masst);
    int n_i = 5;
    int n_f = this->getNf(Q);

    if (n_f < 4) {
        Logger::getInstance()->warn("Scale for alpha_s calculation is below charm mass.");
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

/**
 * @brief Computes the value of the given quark mass at scale Qfin in the MSbar scheme 
 * 
 * @param quark_mass Initial value for the quark mass
 * @param Qinit Scale at which quark_mass is given
 * @param Qfin Scale at which quark_mass should be run
 * @param option_massb Toggle to choose whether to take the running or pole bottom mass
 * @param option_masst Toggle to choose whether to take the running or pole bottom mass
 * @return 
 */
double QCDParameters::running_mass(double quark_mass, double Qinit, double Qfin, std::string option_massb, std::string option_masst) {
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

/**
 * @brief Sets the type of mass (running or pole) to be taken for the calculations
 * @param m_b_type Toggle for b mass
 * @param m_t_type Toggle for t mass
 */
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

/**
 * @brief Runs the given quark mass between two scales with the same number of active flavors in the MSbar scheme
 * 
 * @param mass Initial quark mass
 * @param Q_i Scale of initial quark mass
 * @param Q_f Scale at which quark mass should be run
 * @param nf Number of active flavors at the given scales
 * @return The value of mass(Q_f) in MSbar
 */
double QCDParameters::runMass(double mass, double Q_i, double Q_f, int nf) {
    return mass * this->R(this->runningAlphasCalculation(Q_f, this->m_b_type, this->m_t_type), nf) 
                    / this->R(this->runningAlphasCalculation(Q_i, this->m_b_type, this->m_t_type), nf);
}

/**
 * @brief Computes the number of active flavors at a given energy scale
 * 
 * @param Q Energy scale
 * @return The number of active flavors at scale Q
 */
int QCDParameters::getNf(double Q) {
    auto masses = this->getOrderedMasses();
    for (size_t i = 0; i < masses.size(); ++i) {
        if (1 - Q / masses.at(i) > 1e-4)
            return i;
    }
    return 6;
}

/**
 * @brief Computes the values of the QCD beta function coefficients
 * 
 * @param nf Number of active flavors
 * @return The values of the first three QCD beta function coefficients
 */
std::tuple<double, double, double> QCDParameters::getBetas(int nf) const {
    double b0 = 11 - 2. * nf / 3;
    double b1 = 51 - 19. * nf / 3;
    double b2 = 2857. - 5033. * nf / 9 + 325. * nf * nf / 27;
    return {b0, b1, b2};
}

/**
 * @brief Computes the values of the quark mass anomalous dimension coefficients
 * 
 * @param nf Number of active flavors
 * @return The values of the first three quark mass anomalous dimension coefficients
 */
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

/**
 * @brief Evaluates alpha_s at a given scale given the appropriate value of Lambda_QCD and the corresponding number of active flavors.
 * 
 * @param Q Energy scale 
 * @param Lambda Value Lambda_QCD at the given energy scale
 * @param nf Number of active flavors at the given energy scale
 * @return 
 */
double QCDParameters::alphasRunning(double Q, double Lambda, int nf) const {
    double r = std::pow(Q / Lambda, 2);
    double L = std::log(r);
    double LL = std::log(L);
    auto [b0, b1, b2] = this->getBetas(nf);
    double b02 = b0 * b0;
    double b12 = b1 * b1;
    return 4 * PI * (1 - 2 * b1 * LL / (b02 * L) + 4 * b12 * (std::pow(LL - .5, 2) + b2 * b0 / 8 / b12 - 1.25) / std::pow(b02 * L, 2)) / (b0 * L);
}

/**
 * @brief Computes the value of Lambda_QCD for a given number of active flavors
 * 
 * @param target_alpha Known value of alpha_s at scale Q
 * @param Q Energy scale
 * @param nf Number of active flavors at energy scale Q
 * @return 
 */
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

/**
 * @brief Computes the bottom quark pole mass 
 * 
 * @return The bottom quark pole mass
 */
double QCDParameters::mb_pole() {
	double alphas_mb = runningAlphasCalculation(this->mass_b_b, "running", "pole");	
    this->mass_b_pole = this->mass_b_b * (1. + alphas_mb / PI * (4. / 3. 
            + alphas_mb / PI * ((13.4434 - 1.0414 * 4. + 1.0414 * 4. / 3. * ((this->mass_u + this->mass_d + this->mass_s + this->mass_c) / this->mass_b_b)))));
 	return this->mass_b_pole;
}

/**
 * @brief Computes the charm quark pole mass 
 * 
 * @return The charm quark pole mass
 */
double QCDParameters::mc_pole() {
	double alphas_mc = runningAlphasCalculation(this->mass_c);	
 	return this->mass_c * (1 + alphas_mc / PI * (4. / 3. + alphas_mc / PI * ((13.4434 - 1.0414 * 3. 
            + 1.0414 * 4. / 3. * ((this->mass_u + this->mass_d + this->mass_s) / this->mass_c)))));
}

/**
 * @brief Computes the 1S bottom quark mass 
 * 
 * @return The 1S bottom quark mass
 */
double QCDParameters::mb_1S() {
	double mu = this->mass_b_pole / 2.;
	double alpha = runningAlphasCalculation(mu);
	return this->mass_b_pole * (1 - 2. / 9 * pow(alpha, 2.));
}

/**
 * @brief Computes the top quark running mass at m_top in MSbar
 * 
 * @return The top quark running mass at m_top in MSbar
 */
double QCDParameters::mt_mt() {
	double alpha = runningAlphasCalculation(this->mass_t_pole, "running"); 
    double a = 307. / 32 + PI2 / 3. + PI2 / 9. * log(2) - 1. / 6 * ZETA3 - 71. / 144 * 5;
    this->mass_t_t = this->mass_t_pole / (1 + alpha / PI * (4. / 3 + alpha / PI * a));
    alpha = runningAlphasCalculation(this->mass_t_t, "running", "running");
    this->mass_t_t = this->mass_t_pole / (1 + alpha / PI * (4. / 3 + alpha / PI * a));
	return mass_t_t;
}
