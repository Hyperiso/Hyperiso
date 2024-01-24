#include <vector>
#include "QCDParameters.h"


class Parameters {
public:
    
    double A_b, tan_beta, mu_Q, mass_gluino, mass_b1, mass_b2;
    double inv_alpha_em, M2_Q, mass_t1, mass_t2, A_t, MqL3_Q, MbR_Q, mass_stl, mass_cha1, mass_cha2;
    std::vector<double> yut, yub, mass_neut;
    std::vector<std::vector<double>> sbot_mix, charg_Umix, charg_Vmix, stop_mix, neut_mix;
    std::vector<std::vector<double>> stop_tan_betamix;
    SM sm;
    QCDParameters run;
    double Q {sm.mass_top_pole};

    static Parameters* GetInstance();
    void setScale(double Q); 
private:
    static Parameters* instance;
    Parameters(); // Constructeur pour initialiser les paramètres

    Parameters(const Parameters&) = delete;
    Parameters& operator=(const Parameters&) = delete;
    Parameters(Parameters&&) noexcept = default;
    Parameters& operator=(Parameters&&) noexcept = default;

    
};

struct SM {
    double SM, gp, g2, MSOFT_Q, mass_top_pole, mass_b_pole, mass_b_Q, mass_t_Q;
};
