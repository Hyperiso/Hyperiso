#include "epsilon_calculator.h"
#include "Math.h"


//Ive put MSOFT to the block MSOFT and 1

EpsilonCalculator::EpsilonCalculator() {}


double EpsilonCalculator::epsilon_0() {

    double sw2 = std::pow(std::sin(std::atan((*sm)("GAUGE",1)/ (*sm)("GAUGE",2))), 2);
    double alphas_MSOFT = sm->alpha_s(susy->get_susy_Q());

    

    double term1 = 2.0 / 3.0 * alphas_MSOFT / M_PI * (((*susy)("AD",33) / (*susy)("HMIX",2) - mu_Q) / (*susy)("MASS",1000021) *
               H2((*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)));
    double term2 = -0.5 * (B((*susy)("MASS",1000021), (*susy)("MASS",1000005), susy->get_susy_Q()) + B((*susy)("MASS",1000021), (*susy)("MASS",2000005), susy->get_susy_Q())) / (*susy)("HMIX",2);
    double term3 = 1.0 / (*sm)("SMINPUTS",1) / sw2 / 4.0 / M_PI * (mu_Q * (*susy)("MSOFT",2)) * 
               (((*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) * H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005), mu_Q * mu_Q / (*susy)("MASS",1000005) / (*susy)("MASS",1000005)) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005) / 2.0) +
               ((*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) * H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005), mu_Q * mu_Q / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005) / 2.0));

    Logger *logger = Logger::getInstance();
    logger->info("term1 is " + std::to_string(term3));
    logger->info("sminput is " + std::to_string((*sm)("SMINPUTS",1)));
    return term1 + term2 + term3; 
}


double EpsilonCalculator::epsilon_2() const {

    double sw2 = std::pow(std::sin(std::atan((*sm)("GAUGE",1)/ (*sm)("GAUGE",2))), 2);



    double term1 = (*susy)("YU", 33) * (*susy)("YU", 33) / 16.0 / M_PI / M_PI * 
                   (mu_Q / (*susy)("HMIX",2) - (*susy)("AU", 33)) * 
                   (((*susy)("UMIX",12) * (*susy)("VMIX",12) / (*susy)("MASS",1000024) * 
                     H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000024) / (*susy)("MASS",1000024), (*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000024) / (*susy)("MASS",1000024))) +
                    ((*susy)("UMIX",22) * (*susy)("VMIX",22) / (*susy)("MASS",1000037) * 
                     H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000037) / (*susy)("MASS",1000037), (*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000037) / (*susy)("MASS",1000037))));
    
    double term2 = 1.0 / (*sm)("SMINPUTS",1) / sw2 / 4.0 / M_PI * (mu_Q * (*susy)("MSOFT",2)) * 
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

    double sw2 = std::pow(std::sin(std::atan((*sm)("GAUGE",1)/ (*sm)("GAUGE",2))), 2);
    double alphas_MSOFT = (*sm).QCDRunner.runningAlphasCalculation(susy->get_susy_Q());
    int nb_neut = ((*susy)("MASS", 1000039) == 0.) ? 4 : 5; //mass_neut[5] is gravitino ?

    double epsilonbp = 2.0 / 3.0 * alphas_MSOFT / M_PI * 
                       ((*susy)("AD",33)/ (*susy)("HMIX",2) - mu_Q) / (*susy)("MASS",1000021) * 
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
                     ((*susy)("AU", 33) - mu_Q / (*susy)("HMIX",2)) / (*susy)("MASS",neutralino[ie]) *
                     ((*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * (*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) *
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * (*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) *
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",12) * (*susy)("STOPMIX",12) * (*susy)("SBOTMIX",11) * (*susy)("SBOTMIX",11) *
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",12)* (*susy)("STOPMIX",12) * (*susy)("SBOTMIX",12) * (*susy)("SBOTMIX",12) *
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])));
    }

    epsilonbp += 1.0 / (*sm)("SMINPUTS",1) / sw2 / 4.0 / M_PI * 
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

    double alphas_MSOFT = (*sm).QCDRunner.runningAlphasCalculation(susy->get_susy_Q());
    int nb_neut = ((*susy)("MASS", 1000039) == 0.) ? 4 : 5;

    double epsilon0p = -2.0 / 3.0 * alphas_MSOFT / M_PI * 
                       (mu_Q + (*susy)("AU", 33) / (*susy)("HMIX",2)) / (*susy)("MASS",1000021) *
                       ((*susy)("STOPMIX",11) * (*susy)("STOPMIX",11) * 
                        H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000003) * (*susy)("MASS",1000003) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
                        (*susy)("STOPMIX",12) * (*susy)("STOPMIX",12) * 
                        H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000003) * (*susy)("MASS",1000003) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)));

    for(int ie = 1; ie <= nb_neut; ++ie) {
        epsilon0p += (*susy)("YD", 33) * (*susy)("YD", 33) / 16.0 / M_PI / M_PI * 
                     (*susy)("NMIX", ie*10+4) * (*susy)("NMIX", ie*10+3) * 
                     (mu_Q / (*susy)("HMIX",2)) / (*susy)("MASS",neutralino[ie]) *
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
                   ((*susy)("YD", 33) * (*susy)("YD", 33) * (*susy)("AD",33) / mu_Q * 
                    H2(std::pow((*susy)("MSOFT",43) / mu_Q, 2), std::pow((*susy)("MSOFT", 49) / mu_Q, 2))); //MbR_Q

    // Calcul du deuxième terme en utilisant g2, M2_Q, MqL3_Q, mu_Q
    double term2 = -(*sm)("GAUGE", 2) * (*sm)("GAUGE", 2) * (*susy)("MSOFT",2) / mu_Q * 
                   H2(std::pow((*susy)("MSOFT",43) / mu_Q, 2), std::pow((*susy)("MSOFT",2) / mu_Q, 2)) / 16.0 / M_PI / M_PI;

    return term1 + term2;
}

EpsilonCalculator* EpsilonCalculator::instance = nullptr;