#include "epsilon_calculator.h"
#include "../Math/Math.h"



EpsilonCalculator::EpsilonCalculator() {}


double EpsilonCalculator::epsilon_0() {

    double sw2 = std::pow(std::sin(std::atan((*sm)("Coupling",1)/ (*sm)("Coupling",2))), 2);
    double alphas_MSOFT = (*sm).run.runningAlphasCalculation((*susy).MSOFT_Q);

    

    // Supposons que les valeurs suivantes sont définies dans param
    // Exemple: param.A_b, param.tan_beta, mu_Q, (*susy)("MASS",1000021), etc.

    double term1 = 2.0 / 3.0 * alphas_MSOFT / M_PI * (((*susy)("EXTPAR",12) / (*susy)("EXTPAR",25) - mu_Q) / (*susy)("MASS",1000021) *
               H2((*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)));
    double term2 = -0.5 * (B((*susy)("MASS",1000021), (*susy)("MASS",1000005), (*susy).MSOFT_Q) + B((*susy)("MASS",1000021), (*susy)("MASS",2000005), (*susy).MSOFT_Q)) / (*susy)("EXTPAR",25);
    double term3 = 1.0 / (*sm).inv_alpha_em / sw2 / 4.0 / M_PI * (mu_Q * (*susy)("MSOFT",2)) * 
               (((*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) * H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005), mu_Q * mu_Q / (*susy)("MASS",1000005) / (*susy)("MASS",1000005)) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005) / 2.0) +
               ((*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) * H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005), mu_Q * mu_Q / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005) / 2.0));

    return term1 + term2 + term3;
}


double EpsilonCalculator::epsilon_2() const {

    double sw2 = std::pow(std::sin(std::atan((*sm)("Coupling",1)/ (*sm)("Coupling",2))), 2);

    // Supposons que les valeurs suivantes sont définies dans param
    // Exemple: (*susy)("YU", 33), mu_Q, param.tan_beta, (*susy)("EXTPAR", 11), etc.

    double term1 = (*susy)("YU", 33) * (*susy)("YU", 33) / 16.0 / M_PI / M_PI * 
                   (mu_Q / (*susy)("EXTPAR",25) - (*susy)("EXTPAR", 11)) * 
                   (((*susy)("UMIX",12) * (*susy)("VMIX",12) / (*susy)("MASS",1000024) * 
                     H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000024) / (*susy)("MASS",1000024), (*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000024) / (*susy)("MASS",1000024))) +
                    ((*susy)("UMIX",22) * (*susy)("VMIX",22) / (*susy)("MASS",1000037) * 
                     H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000037) / (*susy)("MASS",1000037), (*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000037) / (*susy)("MASS",1000037))));
    
    double term2 = 1.0 / (*sm).inv_alpha_em / sw2 / 4.0 / M_PI * (mu_Q * (*susy)("MSOFT",2)) * 
                   (((*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * 
                     H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006), mu_Q * mu_Q / (*susy)("MASS",1000006) / (*susy)("MASS",1000006)) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006)) +
                    ((*susy)("STOPMIX",12) * (*susy)("STOPMIX",12)* 
                     H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006), mu_Q * mu_Q / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)));

    return term1 + term2;
}

double EpsilonCalculator::epsilon_b() {
    return epsilon_0() + epsilon_2();
}

// Poursuite de EpsilonCalculator.cpp

// Implémentation de epsilon_bp
double EpsilonCalculator::epsilon_bp() {

    double sw2 = std::pow(std::sin(std::atan((*sm)("Coupling",1)/ (*sm)("Coupling",2))), 2);
    double alphas_MSOFT = (*sm).run.runningAlphasCalculation((*susy).MSOFT_Q);
    int nb_neut = ((*sm).mass_neut[5] == 0.) ? 4 : 5;

    double epsilonbp = 2.0 / 3.0 * alphas_MSOFT / M_PI * 
                       ((*susy)("EXTPAR",12)/ (*susy)("EXTPAR",25) - mu_Q) / (*susy)("MASS",1000021) * 
                       ((*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * (*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11)*
                        H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
                        (*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * (*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) *
                        H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
                        (*susy)("STOPMIX",12) * (*susy)("STOPMIX",12) * (*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) *
                        H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
                        (*susy)("STOPMIX",12) * (*susy)("STOPMIX",12) * (*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) *
                        H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)));

    for(int ie = 1; ie <= nb_neut; ++ie) {
        epsilonbp += (*susy)("YU", 33) * (*susy)("YU", 33) / 16.0 / M_PI / M_PI * 
                     (*susy)("NMIX", ie*10+4) * (*susy)("NMIX", ie*10+3) * 
                     ((*susy)("EXTPAR", 11) - mu_Q / (*susy)("EXTPAR",25)) / (*susy)("MASS",neutralino[ie]) *
                     ((*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * (*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) *
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * (*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) *
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",12) * (*susy)("STOPMIX",12) * (*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) *
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",12)* (*susy)("STOPMIX",12) * (*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) *
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])));
    }

    epsilonbp += 1.0 / (*sm).inv_alpha_em / sw2 / 4.0 / M_PI * 
                 (mu_Q * (*susy)("MSOFT",2)) * 
                 (((*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) *
                   H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006), mu_Q * mu_Q / (*susy)("MASS",1000006) / (*susy)("MASS",1000006)) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006) +
                   (*susy)("STOPMIX",12) * (*susy)("STOPMIX",12) *
                   H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006), mu_Q * mu_Q / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)) / 2.0 +
                 ((*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) *
                  H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005), mu_Q * mu_Q / (*susy)("MASS",1000005) / (*susy)("MASS",1000005)) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005) +
                  (*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) *
                  H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005), mu_Q * mu_Q / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)));

    return epsilonbp;
}

// Poursuite de EpsilonCalculator.cpp

// Implémentation de epsilon_0p
double EpsilonCalculator::epsilon_0p() {

    double alphas_MSOFT = (*sm).run.runningAlphasCalculation((*susy).MSOFT_Q);
    int nb_neut = ((*susy)("MASS", 1000039) == 0.) ? 4 : 5;

    double epsilon0p = -2.0 / 3.0 * alphas_MSOFT / M_PI * 
                       (mu_Q + (*susy)("EXTPAR", 11) / (*susy)("EXTPAR",25)) / (*susy)("MASS",1000021) *
                       ((*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * 
                        H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("EXTPAR",42) * (*susy)("EXTPAR",42) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
                        (*susy)("STOPMIX",12) * (*susy)("STOPMIX",12) * 
                        H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("EXTPAR",42) * (*susy)("EXTPAR",42) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)));

    for(int ie = 1; ie <= nb_neut; ++ie) {
        epsilon0p += (*susy)("YD", 33) * (*susy)("YD", 33) / 16.0 / M_PI / M_PI * 
                     (*susy)("NMIX", ie*10+4) * (*susy)("NMIX", ie*10+3) * 
                     (mu_Q / (*susy)("EXTPAR",25)) / (*susy)("MASS",neutralino[ie]) *
                     ((*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * (*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) * 
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * (*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) * 
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",12) * (*susy)("STOPMIX",12) * (*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) * 
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",12) * (*susy)("STOPMIX",12) * (*susy)("SBOTMIX",12)* (*susy)("SBOTMIX",12) * 
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])));

    }

    return epsilon0p;
}


// Poursuite de EpsilonCalculator.cpp

// Implémentation de epsilon_1p
double EpsilonCalculator::epsilon_1p() const {

    // Calcul du premier terme en utilisant yub[3], A_b, MqL3_Q, MbR_Q, mu_Q
    double term1 = 1.0 / 16.0 / M_PI / M_PI * 
                   ((*susy)("YD", 33) * (*susy)("YD", 33) * (*susy)("EXTPAR",12) / mu_Q * 
                    H2(std::pow((*susy)("EXTPAR",43) / mu_Q, 2), std::pow((*susy)("EXTPAR", 49) / mu_Q, 2))); //MbR_Q

    // Calcul du deuxième terme en utilisant g2, M2_Q, MqL3_Q, mu_Q
    double term2 = -(*sm)("COUPLING", 2) * (*sm)("COUPLING", 2) * (*susy)("MSOFT",2) / mu_Q * 
                   H2(std::pow((*susy)("EXTPAR",43) / mu_Q, 2), std::pow((*susy)("MSOFT",2) / mu_Q, 2)) / 16.0 / M_PI / M_PI;

    return term1 + term2;
}

EpsilonCalculator* EpsilonCalculator::instance = nullptr;