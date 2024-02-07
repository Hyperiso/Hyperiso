#pragma once
#include <vector>
#include "./../Physical_Model/QCDParameters.h"
#include <map>
#include <string>

class Parameters {
public:
    
    double A_b, tan_beta, mu_Q, mass_gluino, mass_b1, mass_b2;
    double inv_alpha_em, M2_Q, mass_t1, mass_t2, A_t, MqL3_Q, MbR_Q, mass_stl, mass_cha1, mass_cha2;
    std::vector<double> yut, yub, mass_neut;
    std::vector<std::vector<double>> sbot_mix, charg_Umix, charg_Vmix, stop_mix, neut_mix;
    std::vector<std::vector<double>> stop_tan_betamix;
    std::vector<std::vector<double>> lambda_u, lambda_d;
    // SM sm;
    
    QCDParameters run;
    // double Q {sm.mass_top_pole};

    static Parameters* GetInstance();
    void setScale(double Q);

    double operator()(std::string block, int pdgCode) {
        if (block == "MASS") {
            return masses[pdgCode];
        }
        if (block =="Coupling") {
            return coupling[pdgCode];
        }
        if (block=="YUKAWA_CH_U") {
            return lambda_u[pdgCode/10][pdgCode%10];
        }
        if (block=="YUKAWA_CH_D") {
            return lambda_d[pdgCode/10][pdgCode%10];
        }
        if (block=="EXTPAR") {
            return extpar[pdgCode];
        }
        return 0;
        
    }
private:
    static Parameters* instance;
    Parameters(); // Constructeur pour initialiser les paramètres

    std::map<int, double> masses;
    std::map<int, double> coupling;
    std::map<int, double> minpar;

    Parameters(const Parameters&) = delete;
    Parameters& operator=(const Parameters&) = delete;
    Parameters(Parameters&&) noexcept = default;
    Parameters& operator=(Parameters&&) noexcept = default;

    
    
};

// struct SM {
//     double SM, gp, g2, MSOFT_Q, mass_top_pole, mass_b_pole, mass_b_Q, mass_t_Q;
// };
