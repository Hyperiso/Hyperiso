#pragma once
#include <vector>
#include "./../Physical_Model/QCDParameters.h"
#include <map>
#include <string>

class Parameters {
public:
    
    double A_b, tan_beta, mu_Q, mass_gluino, mass_b1, mass_b2;
    double inv_alpha_em, M2_Q, mass_t1, mass_t2, A_t, MqL3_Q, MbR_Q, mass_stl, mass_cha1, mass_cha2;
    double MSOFT_Q;
    double inv_alpha_em;
    std::vector<double> yut, yub, mass_neut;
    std::vector<std::vector<double>> sbot_mix, charg_Umix, charg_Vmix, stop_mix, neut_mix;
    std::vector<std::vector<double>> stop_tan_betamix;
    std::vector<std::vector<double>> lambda_u, lambda_d;
    // SM sm;
    
    QCDParameters run;
    // double Q {sm.mass_top_pole};

    static Parameters* GetInstance(int index = 0);
    static Parameters* GetInstance(int index = 0);
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
        if (block=="ALPHA") {
            return alpha[pdgCode];
        }
        if (block=="HMIX") {
            return hmix[pdgCode/10][pdgCode%10];
        }
        if (block=="AMIX") {
            return amix[pdgCode/10][pdgCode%10];
        }
        if (block == "CKM") {
            return ckm[pdgCode/10][pdgCode%10];
        }
        if (block == "STOPMIX") {
            return stopmix[pdgCode/10][pdgCode%10];
        }
        if (block == "SBOTMIX") {
            return stopmix[pdgCode/10][pdgCode%10];
        }
        if (block == "UMIX") {
            return umix[pdgCode/10][pdgCode%10];
        }
        if (block == "VMIX") {
            return vmix[pdgCode/10][pdgCode%10];
        }
        if (block == "MSOFT") {
            return msoft[pdgCode];
        }
        if (block == "YU") {
            return yu[pdgCode/10][pdgCode%10];
        }
        if (block == "YD") {
            return yd[pdgCode/10][pdgCode%10];
        }
        if (block== "NMIX") {
            return nmix[pdgCode/10][pdgCode%10];
        }
        return 0;
        
    }
private:
    static Parameters* instance[2];
    Parameters(); // Constructeur pour initialiser les paramètres

    std::map<int, double> masses;
    std::map<int, double> coupling;
    std::map<int, double> minpar;
    std::map<int, double> extpar;
    std::map<int, double> alpha;
    std::map<int, double> msoft;

    std::vector<std::vector<double>> hmix;
    std::vector<std::vector<double>> amix;
    std::vector<std::vector<double>> stopmix;
    std::vector<std::vector<double>> sbotmix;
    std::vector<std::vector<double>> umix;
    std::vector<std::vector<double>> vmix;
    std::vector<std::vector<double>> nmix;
    std::vector<std::vector<double>> yu;
    std::vector<std::vector<double>> yd;

    std::vector<std::vector<double>> ckm;
    Parameters(const Parameters&) = delete;
    Parameters& operator=(const Parameters&) = delete;
    Parameters(Parameters&&) noexcept = default;
    Parameters& operator=(Parameters&&) noexcept = default;

    
    
};

// struct SM {
//     double SM, gp, g2, MSOFT_Q, mass_top_pole, mass_b_pole, mass_b_Q, mass_t_Q;
// };
