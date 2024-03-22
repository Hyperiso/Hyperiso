#pragma once
#include <cmath>
#include <string>
#include "Math.h"
class QCDParameters {
public:
    // Paramètres à initialiser
    int nf {5};
    double mass_Z;        // MZ
    double alphas_MZ;     // alphas(MZ)
    double Lambda5 {0.};
    double Lambda3, Lambda4, Lambda6;
    double alphasMZ_Lambda3, alphasMZ_Lambda4, alphasMZ_Lambda5, alphasMZ_Lambda6;
    double mass_c;
    double mass_b {4.19999981};
    double mass_t {172.399994};
    double mass_s;
    const double pi {3.14159265358979323846};

    double mass_t_pole;
    double mass_b_pole;
    double mass_b_b;
    double mass_t_t;
    
    QCDParameters() {mass_Z = 91.1699982; alphas_MZ = 0.117200002;}
    QCDParameters(double alpha_Z, double m_Z, double masst_pole, double massb_b, double mass_u, double mass_d, double mass_s, double mass_c);
    QCDParameters& operator=(const QCDParameters& other) {
        if (this != &other) {
            this->mass_t_pole = other.mass_t_pole;
            this->mass_b_pole = other.mass_b_pole;
            this->mass_b_b = other.mass_b_b;
            this->mass_t_t = other.mass_t_t;
        }
        return *this;
    }

    double alphasRunning(double Q, double Lambda, int nf) const;
    double matchLambda(double alpha_running, double Q, int nf);
    double runningAlphasCalculation(double Q, std::string option_massb = "pole", std::string option_masst = "pole");
    double running_mass(double quark_mass, double Qinit, double Qfin,  std::string option_massb = "pole", std::string option_masst = "pole");

    double mb_pole(double mass_b, double mass_u, double mass_d, double mass_s, double mass_c);
    double mc_pole(double mass_u, double mass_d, double mass_s, double mass_c);
    double mb_1S(double mb_pole);
    double mt_mt(double mt_pole);

    double get_mt_mt() {return this->mass_t_t;}
};
