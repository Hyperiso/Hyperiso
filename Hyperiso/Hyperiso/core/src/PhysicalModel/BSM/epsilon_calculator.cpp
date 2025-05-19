#include "epsilon_calculator.h"
#include "Math.h"


//Ive put MSOFT to the block MSOFT and 1

// EpsilonCalculator::EpsilonCalculator() {}


// double EpsilonCalculator::epsilon_0() {

//     double sw2 = std::pow(std::sin(std::atan((*sm)("GAUGE",1)/ (*sm)("GAUGE",2))), 2);
//     double alphas_MSOFT = QCDHelper::alpha_s((*susy)("HMIX", 0));

//     double factor = 2.0 / 3.0 * alphas_MSOFT / M_PI;

//     double term1 =  ((*susy)("AD",22) / (*susy)("HMIX",2) - mu_Q) / (*susy)("MASS",1000021) *
//                H2((*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021));
//     double term2 = -0.5 * (B((*susy)("MASS",1000021), (*susy)("MASS",1000005), (*susy)("HMIX", 0)) + B((*susy)("MASS",1000021), (*susy)("MASS",2000005), (*susy)("HMIX", 0))) / (*susy)("HMIX",2);
//     double term3 = 1.0 / (*sm)("SMINPUTS",1) / sw2 / 4.0 / M_PI * (mu_Q * (*susy)("MSOFT",2)) * 
//                ((*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) * H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005), mu_Q * mu_Q / (*susy)("MASS",1000005) / (*susy)("MASS",1000005)) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005) / 2.0 +
//                (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) * H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005), mu_Q * mu_Q / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005) / 2.0);

//     LOG_DEBUG("AD 22 " + std::to_string((*susy)("AD",22)));
//     LOG_DEBUG("term1 in epsilon_0 is " + std::to_string(term1));
//     LOG_DEBUG("term2 in epsilon_0 is " + std::to_string(term2));
//     LOG_DEBUG("term3 in epsilon_0 is " + std::to_string(term3));
//     return factor * (term1 + term2) + term3; 
// }


// double EpsilonCalculator::epsilon_2() const {

//     double sw2 = std::pow(std::sin(std::atan((*sm)("GAUGE",1)/ (*sm)("GAUGE",2))), 2);



//     double term1 = (*susy)("YU", 22) * (*susy)("YU", 22) / 16.0 / M_PI / M_PI * 
//                    (mu_Q / (*susy)("HMIX",2) - (*susy)("AU", 22)) * 
//                    (((*susy)("UMIX",01) * (*susy)("VMIX",01) / (*susy)("MASS",1000024) * 
//                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000024) / (*susy)("MASS",1000024), (*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000024) / (*susy)("MASS",1000024))) +
//                     ((*susy)("UMIX",11) * (*susy)("VMIX",11) / (*susy)("MASS",1000037) * 
//                      H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000037) / (*susy)("MASS",1000037), (*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000037) / (*susy)("MASS",1000037))));
    
//     double term2 = 1.0 / (*sm)("SMINPUTS",1) / sw2 / 4.0 / M_PI * (mu_Q * (*susy)("MSOFT",2)) * 
//                    (((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * 
//                      H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006), mu_Q * mu_Q / (*susy)("MASS",1000006) / (*susy)("MASS",1000006)) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006)) +
//                     ((*susy)("STOPMIX",01) * (*susy)("STOPMIX",01)* 
//                      H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006), mu_Q * mu_Q / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)));

//     return term1 + term2;
// }

// double EpsilonCalculator::epsilon_b() {

//     LOG_DEBUG("epsilon 0 : " + std::to_string(epsilon_0()));
//     LOG_DEBUG("epsilon 2 : " + std::to_string(epsilon_2()));
//     return epsilon_0() + epsilon_2();
// }

// // Poursuite de EpsilonCalculator.cpp

// // Implémentation de epsilon_bp
// double EpsilonCalculator::epsilon_bp() {

//     double sw2 = std::pow(std::sin(std::atan((*sm)("GAUGE",1)/ (*sm)("GAUGE",2))), 2);
//     double alphas_MSOFT = QCDHelper::alpha_s((*susy)("HMIX", 0));
//     int nb_neut = ((*susy)("MASS", 1000039) == 0.) ? 4 : 5; //mass_neut[5] is gravitino ?

    
//     double epsilonbp = 2.0 / 3.0 * alphas_MSOFT / M_PI * 
//                        ((*susy)("AD",22)/ (*susy)("HMIX",2) - mu_Q) / (*susy)("MASS",1000021) * 
//                        ((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00)*
//                         H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
//                         (*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
//                         H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
//                         (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) *
//                         H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
//                         (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
//                         H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)));
//     LOG_INFO("AD", (*susy)("AD",22));
//     LOG_INFO("stopmix00", (*susy)("STOPMIX",00));
//     LOG_INFO("stopmix01", (*susy)("STOPMIX",01));
//     LOG_INFO("sbopmix00", (*susy)("SBOTMIX",00));
//     LOG_INFO("sbopmix01", (*susy)("SBOTMIX",01));
//     LOG_INFO("t1", (*susy)("MASS",1000006) );
//     LOG_INFO("t2", (*susy)("MASS",2000006) );
//     for(int ie = 0; ie < nb_neut; ++ie) {
//         epsilonbp += (*susy)("YU", 22) * (*susy)("YU", 22) / 16.0 / M_PI / M_PI * 
//                      (*susy)("NMIX", ie*10+3) * (*susy)("NMIX", ie*10+2) * 
//                      ((*susy)("AU", 22) - mu_Q / (*susy)("HMIX",2)) / (*susy)("MASS",neutralino[ie]) *
//                      ((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) *
//                       H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
//                       (*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
//                       H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
//                       (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) *
//                       H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
//                       (*susy)("STOPMIX",01)* (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
//                       H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])));
//     }

//     epsilonbp += 1.0 / (*sm)("SMINPUTS",1) / sw2 / 4.0 / M_PI * 
//                  (mu_Q * (*susy)("MSOFT",2)) * 
//                  (((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) *
//                    H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006), mu_Q * mu_Q / (*susy)("MASS",1000006) / (*susy)("MASS",1000006)) / (*susy)("MASS",1000006) / (*susy)("MASS",1000006) +
//                    (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) *
//                    H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006), mu_Q * mu_Q / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)) / (*susy)("MASS",2000006) / (*susy)("MASS",2000006)) / 2.0 +
//                  ((*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) *
//                   H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005), mu_Q * mu_Q / (*susy)("MASS",1000005) / (*susy)("MASS",1000005)) / (*susy)("MASS",1000005) / (*susy)("MASS",1000005) +
//                   (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) *
//                   H2((*susy)("MSOFT",2) * (*susy)("MSOFT",2) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005), mu_Q * mu_Q / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)) / (*susy)("MASS",2000005) / (*susy)("MASS",2000005)));

//     return epsilonbp;
// }


// double EpsilonCalculator::epsilon_0p() {

//     double alphas_MSOFT = QCDHelper::alpha_s((*susy)("HMIX", 0));
//     int nb_neut = ((*susy)("MASS", 1000039) == 0.) ? 4 : 5;

//     double epsilon0p = -2.0 / 3.0 * alphas_MSOFT / M_PI * 
//                        (mu_Q + (*susy)("AU", 22) / (*susy)("HMIX",2)) / (*susy)("MASS",1000021) *
//                        ((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * 
//                         H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000003) * (*susy)("MASS",1000003) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)) +
//                         (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * 
//                         H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021), (*susy)("MASS",1000003) * (*susy)("MASS",1000003) / (*susy)("MASS",1000021) / (*susy)("MASS",1000021)));

//     for(int ie = 0; ie < nb_neut; ++ie) {
//         epsilon0p += (*susy)("YD", 22) * (*susy)("YD", 22) / 16.0 / M_PI / M_PI * 
//                      (*susy)("NMIX", ie*10+3) * (*susy)("NMIX", ie*10+2) * 
//                      (mu_Q / (*susy)("HMIX",2)) / (*susy)("MASS",neutralino[ie]) *
//                      ((*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) * 
//                       H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
//                       (*susy)("STOPMIX",00) * (*susy)("STOPMIX",00) * (*susy)("SBOTMIX",01) * (*susy)("SBOTMIX",01) * 
//                       H2((*susy)("MASS",1000006) * (*susy)("MASS",1000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
//                       (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",00) * (*susy)("SBOTMIX",00) * 
//                       H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",2000005) * (*susy)("MASS",2000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])) +
//                       (*susy)("STOPMIX",01) * (*susy)("STOPMIX",01) * (*susy)("SBOTMIX",01)* (*susy)("SBOTMIX",01) * 
//                       H2((*susy)("MASS",2000006) * (*susy)("MASS",2000006) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie]), (*susy)("MASS",1000005) * (*susy)("MASS",1000005) / (*susy)("MASS",neutralino[ie]) / (*susy)("MASS",neutralino[ie])));

//     }

//     return epsilon0p;
// }


// // Poursuite de EpsilonCalculator.cpp

// // Implémentation de epsilon_1p
// double EpsilonCalculator::epsilon_1p() const {

//     // Calcul du premier terme en utilisant yub[3], A_b, MqL3_Q, MbR_Q, mu_Q
//     double term1 = 1.0 / 16.0 / M_PI / M_PI * 
//                    ((*susy)("YD", 22) * (*susy)("YD", 22) * (*susy)("AD",22) / mu_Q * 
//                     H2(std::pow((*susy)("MSOFT",43) / mu_Q, 2), std::pow((*susy)("MSOFT", 49) / mu_Q, 2))); //MbR_Q

//     // Calcul du deuxième terme en utilisant g2, M2_Q, MqL3_Q, mu_Q
//     double term2 = -(*sm)("GAUGE", 2) * (*sm)("GAUGE", 2) * (*susy)("MSOFT",2) / mu_Q * 
//                    H2(std::pow((*susy)("MSOFT",43) / mu_Q, 2), std::pow((*susy)("MSOFT",2) / mu_Q, 2)) / 16.0 / M_PI / M_PI;

//     return term1 + term2;
// }


void EpsilonCalculator::init() {
    if (EpsilonCalculator::initialized) {
		return;
	}

    init_epsilon_block();
}

void EpsilonCalculator::init_epsilon_block() {
    std::unordered_map<ParameterType, std::vector<std::string>> src = {
        {ParameterType::SM, {"MASS", "SMINPUTS"}}, 
        {ParameterType::BSM, {"MASS", "GAUGE", "HMIX", "MSOFT", "AD", "AU", "YD", "YU", "SBOTMIX", "STOPMIX", "UMIX", "VMIX", "NMIX", "ALPHA"}},
        {ParameterType::WILSON, {"WPARAM_SI_SM"}}
    };

    auto func = [] (const std::unordered_map<std::string, std::shared_ptr<Block>>& src, std::shared_ptr<DependentBlock> dep_block) {

        src.at("ALPHA")->retrieve({})->get_val();

        double g2 = src.at("GAUGE")->retrieve(2)->get_val();
        double alpha_em = src.at("SMINPUTS")->retrieve(1)->get_val();

        double m_ds = src.at("MASS")->retrieve(1000001)->get_val();
        double m_us = src.at("MASS")->retrieve(1000002)->get_val();
        double m_ss = src.at("MASS")->retrieve(1000003)->get_val();
        double m_cs = src.at("MASS")->retrieve(1000004)->get_val();
        double m_bs = src.at("MASS")->retrieve(1000005)->get_val();
        double m_ts = src.at("MASS")->retrieve(1000006)->get_val();

        double m_d2s = src.at("MASS")->retrieve(2000001)->get_val();
        double m_u2s = src.at("MASS")->retrieve(2000002)->get_val();
        double m_s2s = src.at("MASS")->retrieve(2000003)->get_val();
        double m_c2s = src.at("MASS")->retrieve(2000004)->get_val();
        double m_b2s = src.at("MASS")->retrieve(2000005)->get_val();
        double m_t2s = src.at("MASS")->retrieve(2000006)->get_val();

        double m_gluino = src.at("MASS")->retrieve(1000021)->get_val();
        double m_c1 = src.at("MASS")->retrieve(1000024)->get_val();
        double m_c2 = src.at("MASS")->retrieve(1000037)->get_val();
        double mu_Q = src.at("HMIX")->retrieve(1)->get_val();

        double mqL3 = src.at("MSOFT")->retrieve(43)->get_val();
        double mbR = src.at("MSOFT")->retrieve(49)->get_val();

        double m_G = src.at("MASS")->retrieve(1000039)->get_val(); //TODO

        std::map<int,int> neutralino = {{0, 1000022},{1, 1000023},{2, 1000025},{3, 1000035}};


        std::vector<double> m_neutralino = {src.at("MASS")->retrieve(neutralino[0])->get_val(), src.at("MASS")->retrieve(neutralino[1])->get_val(), src.at("MASS")->retrieve(neutralino[2])->get_val(), src.at("MASS")->retrieve(neutralino[3])->get_val()};
        double ad_22 = src.at("AD")->retrieve({3, 3})->get_val();
        double au_22 = src.at("AU")->retrieve({3, 3})->get_val();

        double yu_22 = src.at("YU")->retrieve({3, 3})->get_val();
        double yd_22 = src.at("YD")->retrieve({3, 3})->get_val();

        double sbot_mix_00 = src.at("SBOTMIX")->retrieve({0+1, 0+1})->get_val();
        double sbot_mix_01 = src.at("SBOTMIX")->retrieve({0+1, 1+1})->get_val();

        double stop_mix_00 = src.at("STOPMIX")->retrieve({0+1, 0+1})->get_val();
        double stop_mix_01 = src.at("STOPMIX")->retrieve({0+1, 1+1})->get_val();

        double umix_01 = src.at("UMIX")->retrieve({0+1, 1+1})->get_val();
        double umix_11 = src.at("UMIX")->retrieve({1+1, 1+1})->get_val();

        double vmix_01 = src.at("VMIX")->retrieve({0+1, 1+1})->get_val();
        double vmix_11 = src.at("VMIX")->retrieve({1+1, 1+1})->get_val();
        
        //0
        double sw2 = src.at("WPARAM_SI_SM")->retrieve(4)->get_val();

        // double alphas_MSOFT = QCDHelper::alpha_s(src.at("HMIX")->retrieve(0)->get_val()); // SUSY Breaking scale
        double alphas_MSOFT = QCDHelper::alpha_s(2.448e3);

        double tan_beta = src.at("HMIX")->retrieve(2)->get_val();

        double factor = 2.0 / 3.0 * alphas_MSOFT / M_PI;

        double M_2 = src.at("MSOFT")->retrieve(2)->get_val();


        double term1 =  (ad_22 / tan_beta - mu_Q) / m_gluino *
                H2(m_bs * m_bs / m_gluino / m_gluino, m_b2s * m_b2s / m_gluino / m_gluino);
        double term2 = -0.5 * (B(m_gluino, m_bs, alphas_MSOFT) + B(m_gluino, m_b2s, alphas_MSOFT)) / tan_beta;
        double term3 = 1.0 / alpha_em / sw2 / 4.0 / M_PI * (mu_Q * M_2) * 
                (sbot_mix_00 * sbot_mix_00 * H2(M_2 * M_2 / m_bs / m_bs, mu_Q * mu_Q / m_bs / m_bs) / m_bs / m_bs / 2.0 +
                sbot_mix_01 * sbot_mix_01 * H2(M_2 * M_2 / m_b2s / m_b2s, mu_Q * mu_Q / m_b2s / m_b2s) / m_b2s / m_b2s / 2.0);


        double epsilon_0 = factor * (term1 + term2) + term3;


        term1 = yu_22 * yu_22 / 16.0 / M_PI / M_PI * 
                    (mu_Q / tan_beta - au_22) * 
                    ((umix_01 * vmix_01 / m_c1 * 
                        H2(m_ts * m_ts / m_c1 / m_c1, m_t2s * m_t2s / m_c1 / m_c1)) +
                        (umix_11 * vmix_11 / m_c2 * 
                        H2(m_ts * m_ts / m_c2 / m_c2, m_t2s * m_t2s / m_c2 / m_c2)));
        
        term2 = 1.0 / alpha_em / sw2 / 4.0 / M_PI * (mu_Q * M_2) * 
                    ((stop_mix_00 * stop_mix_00 * 
                        H2(M_2 * M_2 / m_ts / m_ts, mu_Q * mu_Q / m_ts / m_ts) / m_ts / m_ts) +
                        (stop_mix_01 * stop_mix_01* 
                        H2(M_2 * M_2 / m_t2s / m_t2s, mu_Q * mu_Q / m_t2s / m_t2s) / m_t2s / m_t2s));

        double epsilon_2 = term1 + term2;


        //b
        double epsilon_b = epsilon_0 + epsilon_2;

        //bp
        int nb_neut = (m_G == 0.) ? 4 : 5; //mass_neut[5] is gravitino ?

        
        double epsilonbp = 2.0 / 3.0 * alphas_MSOFT / M_PI * 
                        (ad_22/ tan_beta - mu_Q) / m_gluino * 
                        (stop_mix_00 * stop_mix_00 * sbot_mix_00 * sbot_mix_00*
                            H2(m_ts * m_ts / m_gluino / m_gluino, m_b2s * m_b2s / m_gluino / m_gluino) +
                            stop_mix_00 * stop_mix_00 * sbot_mix_01 * sbot_mix_01 *
                            H2(m_ts * m_ts / m_gluino / m_gluino, m_bs * m_bs / m_gluino / m_gluino) +
                            stop_mix_01 * stop_mix_01 * sbot_mix_00 * sbot_mix_00 *
                            H2(m_t2s * m_t2s / m_gluino / m_gluino, m_b2s * m_b2s / m_gluino / m_gluino) +
                            stop_mix_01 * stop_mix_01 * sbot_mix_01 * sbot_mix_01 *
                            H2(m_t2s * m_t2s / m_gluino / m_gluino, m_bs * m_bs / m_gluino / m_gluino));

        for(int ie = 0; ie < nb_neut; ++ie) {
            epsilonbp += yu_22 * yu_22 / 16.0 / M_PI / M_PI * 
                        src.at("NMIX")->retrieve({ie + 1, 3 + 1})->get_val() * src.at("NMIX")->retrieve({ie + 1, 2 + 1})->get_val() * 
                        (au_22 - mu_Q / tan_beta) / m_neutralino[ie] *
                        (stop_mix_00 * stop_mix_00 * sbot_mix_00 * sbot_mix_00 *
                        H2(m_t2s * m_t2s / m_neutralino[ie] / m_neutralino[ie], m_bs * m_bs / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_00 * stop_mix_00 * sbot_mix_01 * sbot_mix_01 *
                        H2(m_t2s * m_t2s / m_neutralino[ie] / m_neutralino[ie], m_b2s * m_b2s / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_01 * stop_mix_01 * sbot_mix_00 * sbot_mix_00 *
                        H2(m_ts * m_ts / m_neutralino[ie] / m_neutralino[ie], m_bs * m_bs / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_01* stop_mix_01 * sbot_mix_01 * sbot_mix_01 *
                        H2(m_ts * m_ts / m_neutralino[ie] / m_neutralino[ie], m_b2s * m_b2s / m_neutralino[ie] / m_neutralino[ie]));
        }

        epsilonbp += 1.0 / alpha_em / sw2 / 4.0 / M_PI * 
                    (mu_Q * M_2) * 
                    ((stop_mix_00 * stop_mix_00 *
                    H2(M_2 * M_2 / m_ts / m_ts, mu_Q * mu_Q / m_ts / m_ts) / m_ts / m_ts +
                    stop_mix_01 * stop_mix_01 *
                    H2(M_2 * M_2 / m_t2s / m_t2s, mu_Q * mu_Q / m_t2s / m_t2s) / m_t2s / m_t2s) / 2.0 +
                    (sbot_mix_00 * sbot_mix_00 *
                    H2(M_2 * M_2 / m_bs / m_bs, mu_Q * mu_Q / m_bs / m_bs) / m_bs / m_bs +
                    sbot_mix_01 * sbot_mix_01 *
                    H2(M_2 * M_2 / m_b2s / m_b2s, mu_Q * mu_Q / m_b2s / m_b2s) / m_b2s / m_b2s));

        // return epsilonbp;
        
        //0p

        double epsilon0p = -2.0 / 3.0 * alphas_MSOFT / M_PI * 
                        (mu_Q + au_22 / tan_beta) / m_gluino *
                        (stop_mix_00 * stop_mix_00 * 
                            H2(m_t2s * m_t2s / m_gluino / m_gluino, m_ss * m_ss / m_gluino / m_gluino) +
                            stop_mix_01 * stop_mix_01 * 
                            H2(m_ts * m_ts / m_gluino / m_gluino, m_ss * m_ss / m_gluino / m_gluino));

        for(int ie = 0; ie < nb_neut; ++ie) {
            epsilon0p += yd_22 * yd_22 / 16.0 / M_PI / M_PI * 
                        src.at("NMIX")->retrieve({ie+1, 3+1})->get_val() * src.at("NMIX")->retrieve({ie+1, 2+1})->get_val() * 
                        (mu_Q / tan_beta) / m_neutralino[ie] *
                        (stop_mix_00 * stop_mix_00 * sbot_mix_00 * sbot_mix_00 * 
                        H2(m_ts * m_ts / m_neutralino[ie] / m_neutralino[ie], m_b2s * m_b2s / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_00 * stop_mix_00 * sbot_mix_01 * sbot_mix_01 * 
                        H2(m_ts * m_ts / m_neutralino[ie] / m_neutralino[ie], m_bs * m_bs / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_01 * stop_mix_01 * sbot_mix_00 * sbot_mix_00 * 
                        H2(m_t2s * m_t2s / m_neutralino[ie] / m_neutralino[ie], m_b2s * m_b2s / m_neutralino[ie] / m_neutralino[ie]) +
                        stop_mix_01 * stop_mix_01 * sbot_mix_01* sbot_mix_01 * 
                        H2(m_t2s * m_t2s / m_neutralino[ie] / m_neutralino[ie], m_bs * m_bs / m_neutralino[ie] / m_neutralino[ie]));

        }

        // return epsilon0p;

        //1p
        // Calcul du premier terme en utilisant yub[3], A_b, MqL3_Q, MbR_Q, mu_Q
        term1 = 1.0 / 16.0 / M_PI / M_PI * 
        (yd_22 * yd_22 * ad_22 / mu_Q * 
        H2(std::pow(mqL3 / mu_Q, 2), std::pow(mbR / mu_Q, 2))); //MbR_Q

        // Calcul du deuxième terme en utilisant g2, M2_Q, MqL3_Q, mu_Q
        term2 = -g2 * g2 * M_2 / mu_Q * 
            H2(std::pow(mqL3 / mu_Q, 2), std::pow(M_2 / mu_Q, 2)) / 16.0 / M_PI / M_PI;

        double epsilon_1p = term1 + term2;
        
        double epsfac=pow((1.+epsilon_b*tan_beta),2.);

        dep_block->store_or_assign({0,1}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", {0,1}}, epsilon_0, 0., 0.));
		dep_block->store_or_assign({0,2}, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", {0,2}}, epsilon0p, 0., 0.));
		dep_block->store_or_assign(1, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 1}, epsilon_1p, 0., 0.));
		dep_block->store_or_assign(2, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 2}, epsilon_2, 0., 0.));
		dep_block->store_or_assign(3, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 3}, epsilon_b, 0., 0.));
        dep_block->store_or_assign(4, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 4}, epsilonbp, 0., 0.));
        dep_block->store_or_assign(5, std::make_shared<Parameter>(ParamId{ParameterType::WILSON, "EPSILON_SUSY", 5}, epsfac, 0., 0.)); 

    };

    EpsilonCalculator::composer.compose_block("EPSILON_SUSY", src, func);
}