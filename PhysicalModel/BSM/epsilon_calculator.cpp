#include "epsilon_calculator.h"
#include "Math.h"


//Ive put MSOFT to the block MSOFT and 1

EpsilonCalculator::EpsilonCalculator() {}


double EpsilonCalculator::epsilon_0() {

    double sw2 = std::pow(std::sin(std::atan((*sm)("GAUGE",1)/ (*sm)("GAUGE",2))), 2);
    double alphas_MSOFT = sm->alpha_s(susy->get_susy_Q());

    double factor = 2.0 / 3.0 * alphas_MSOFT / M_PI;

    double term1 =  ((*susy)("AD",22) / (*susy)("HMIX",2) - mu_Q) / (*susy)("MASS",1000021) *
               H2((*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021));
    double term2 = -0.5 * (B((*susy)("MASS",1000021), (*susy)("MASS",1000005), susy->get_susy_Q()) + B((*susy)("MASS",1000021), (*susy)("MASS",2000005), susy->get_susy_Q())) / (*susy)("HMIX",2);
    double term3 = 1.0 / (*sm)("SMINPUTS",1) / sw2 / 4.0 / M_PI * (mu_Q * (*susy)("MSOFT",2)) * 
               ((*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) * H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005), mu_Q * mu_Q / (*susy)("MASS",1000005) / (*susy)("MASS",1000005)) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005) / 2.0 +
               (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) * H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005), mu_Q * mu_Q / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005) / 2.0);

    Logger *logger = Logger::getInstance();
    LOG_DEBUG("AD 22 " + std::to_string((*susy)("AD",22)));
    LOG_DEBUG("term1 in epsilon_0 is " + std::to_string(term1));
    LOG_DEBUG("term2 in epsilon_0 is " + std::to_string(term2));
    LOG_DEBUG("term3 in epsilon_0 is " + std::to_string(term3));
    return factor * (term1 + term2) + term3; 
}


double EpsilonCalculator::epsilon_2() const {

    double sw2 = std::pow(std::sin(std::atan((*sm)("GAUGE",1)/ (*sm)("GAUGE",2))), 2);



    double term1 = (*susy)("YU", 22) * (*susy)("YU", 22) / 16.0 / M_PI / M_PI * 
                   (mu_Q / (*susy)("HMIX",2) - (*susy)("AU", 22)) * 
                   (((*susy)("UMIX",01) * (*susy)("VMIX",01) / (*susy)("MASS",1000024) * 
                     H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000024) / (*susy)("MASS",1000024), (*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000024) / (*susy)("MASS",1000024))) +
                    ((*susy)("UMIX",11) * (*susy)("VMIX",11) / (*susy)("MASS",1000037) * 
                     H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000037) / (*susy)("MASS",1000037), (*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000037) / (*susy)("MASS",1000037))));
    
    double term2 = 1.0 / (*sm)("SMINPUTS",1) / sw2 / 4.0 / M_PI * (mu_Q * (*susy)("MSOFT",2)) * 
                   (((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * 
                     H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006), mu_Q * mu_Q / (*susy)("MASS",1000006) / (*susy)("MASS",1000006)) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006)) +
                    ((*susy)("STOPMIX",01) * (*susy)("STOPMIX",01)* 
                     H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006), mu_Q * mu_Q / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)));

    return term1 + term2;
}

double EpsilonCalculator::epsilon_b() {

    LOG_DEBUG("epsilon 0 : " + std::to_string(epsilon_0()));
    LOG_DEBUG("epsilon 2 : " + std::to_string(epsilon_2()));
    return epsilon_0() + epsilon_2();
}

// Poursuite de EpsilonCalculator.cpp

// Implémentation de epsilon_bp
double EpsilonCalculator::epsilon_bp() {

    double sw2 = std::pow(std::sin(std::atan((*sm)("GAUGE",1)/ (*sm)("GAUGE",2))), 2);
    double alphas_MSOFT = (*sm).QCDRunner.runningAlphasCalculation(susy->get_susy_Q());
    int nb_neut = ((*susy)("MASS", 1000039) == 0.) ? 4 : 5; //mass_neut[5] is gravitino ?


    double epsilonbp = 2.0 / 3.0 * alphas_MSOFT / M_PI * 
                       ((*susy)("AD",22)/ (*susy)("HMIX",2) - mu_Q) / (*susy)("MASS",1000021) * 
                       ((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00)*
                        H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
                        (*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
                        H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
                        (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) *
                        H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
                        (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
                        H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)));

    for(int ie = 0; ie < nb_neut; ++ie) {
        epsilonbp += (*susy)("YU", 22) * (*susy)("YU", 22) / 16.0 / M_PI / M_PI * 
                     (*susy)("NMIX", ie*10+3) * (*susy)("NMIX", ie*10+2) * 
                     ((*susy)("AU", 22) - mu_Q / (*susy)("HMIX",2)) / (*susy)("MASS",neutralino[ie]) *
                     ((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) *
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) *
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",01)* (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])));
    }

    epsilonbp += 1.0 / (*sm)("SMINPUTS",1) / sw2 / 4.0 / M_PI * 
                 (mu_Q * (*susy)("MSOFT",2)) * 
                 (((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) *
                   H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006), mu_Q * mu_Q / (*susy)("MASS",1000006) / (*susy)("MASS",1000006)) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006) +
                   (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) *
                   H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006), mu_Q * mu_Q / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)) / 2.0 +
                 ((*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) *
                  H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005), mu_Q * mu_Q / (*susy)("MASS",1000005) / (*susy)("MASS",1000005)) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005) +
                  (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
                  H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005), mu_Q * mu_Q / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)));

    return epsilonbp;
}


double EpsilonCalculator::epsilon_0p() {

    double alphas_MSOFT = (*sm).QCDRunner.runningAlphasCalculation(susy->get_susy_Q());
    int nb_neut = ((*susy)("MASS", 1000039) == 0.) ? 4 : 5;

    double epsilon0p = -2.0 / 3.0 * alphas_MSOFT / M_PI * 
                       (mu_Q + (*susy)("AU", 22) / (*susy)("HMIX",2)) / (*susy)("MASS",1000021) *
                       ((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * 
                        H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000003) * (*susy)("MASS",1000003) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
                        (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * 
                        H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000003) * (*susy)("MASS",1000003) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)));

    for(int ie = 0; ie < nb_neut; ++ie) {
        epsilon0p += (*susy)("YD", 22) * (*susy)("YD", 22) / 16.0 / M_PI / M_PI * 
                     (*susy)("NMIX", ie*10+3) * (*susy)("NMIX", ie*10+2) * 
                     (mu_Q / (*susy)("HMIX",2)) / (*susy)("MASS",neutralino[ie]) *
                     ((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) * 
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) * 
                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) * 
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
                      (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",01)* (*susy)("SBOTMIX",01) * 
                      H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])));

    }

    return epsilon0p;
}


// Poursuite de EpsilonCalculator.cpp

// Implémentation de epsilon_1p
double EpsilonCalculator::epsilon_1p() const {

    // Calcul du premier terme en utilisant yub[3], A_b, MqL3_Q, MbR_Q, mu_Q
    double term1 = 1.0 / 16.0 / M_PI / M_PI * 
                   ((*susy)("YD", 22) * (*susy)("YD", 22) * (*susy)("AD",22) / mu_Q * 
                    H2(std::pow((*susy)("MSOFT",43) / mu_Q, 2), std::pow((*susy)("MSOFT", 49) / mu_Q, 2))); //MbR_Q

    // Calcul du deuxième terme en utilisant g2, M2_Q, MqL3_Q, mu_Q
    double term2 = -(*sm)("GAUGE", 2) * (*sm)("GAUGE", 2) * (*susy)("MSOFT",2) / mu_Q * 
                   H2(std::pow((*susy)("MSOFT",43) / mu_Q, 2), std::pow((*susy)("MSOFT",2) / mu_Q, 2)) / 16.0 / M_PI / M_PI;

    return term1 + term2;
}

EpsilonCalculator* EpsilonCalculator::instance = nullptr;