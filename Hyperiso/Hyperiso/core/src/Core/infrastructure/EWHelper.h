#ifndef __EWHELPER_H__
#define __EWHELPER_H__

#include <array>
#include <string>
#include "Parameters.h"
#include "Math.h"
#include "QCDHelper.h"

// All calculations are taken from hep-ph/0011135

class EWHelper {
private:
    static double I_s(int k, double mu_low, double mu_high);
    static double delta_lept(double mu_low, double mu_high, double alpha);
    static double delta_part(double mu_low, double mu_high, int n_f, double alpha);
    static double delta_had(double mu_low, double mu_high, int n_f, double alpha);
    static double pi_heavy(double mu, double e_q, double m_q_pole, double n_l, double alpha);
    static double pi_light_heavy(double mu, double m_q_pole, double n_l, double alpha);

    static double alpha_5_mu_b(double inv_alpha_m_Z, double m_Z, double mu_b);
    static double alpha_4_mu_b(double inv_alpha_m_Z, double m_Z, double m_b_pole, double mu_b);
    static double alpha_4_mu_c(double inv_alpha_m_Z, double m_Z, double m_b_pole, double mu_b, double mu_c);

public:
    static void Init();

    static double alpha_em(double mu);
};

#endif // __EWHELPER_H__
